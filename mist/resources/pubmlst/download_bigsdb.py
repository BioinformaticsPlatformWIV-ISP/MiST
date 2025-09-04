#!/usr/bin/env python3
# Script to download authenticated resources from PubMLST and BIGSdb Pasteur
# via their REST interfaces.
# Written by Keith Jolley
# Copyright (c) 2024, University of Oxford
# E-mail: keith.jolley@biology.ox.ac.uk
#
# BIGSdb_downloader is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BIGSdb_downloader is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# Version 20241030
import argparse
import os
import stat
import re
import configparser
import sys
import json
from urllib.parse import parse_qs
from pathlib import Path
from rauth import OAuth1Service, OAuth1Session

BASE_WEB = {
    "PubMLST": "https://pubmlst.org/bigsdb",
    "Pasteur": "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl",
}
BASE_API = {
    "PubMLST": "https://rest.pubmlst.org",
    "Pasteur": "https://bigsdb.pasteur.fr/api",
}


parser = argparse.ArgumentParser()

parser.add_argument(
    "--cron",
    action="store_true",
    help="Script is being run as a CRON job or non-interactively.",
)
parser.add_argument(
    "--db", required=False, help="Database config - only needed for setup."
)
parser.add_argument(
    "--json_body",
    required=False,
    help="JSON body to be included in a POST call. If this is longer than "
    "the command line limit (probably about 128kb) then you will need to "
    "save the JSON payload to a file and use --json_body_file",
)
parser.add_argument(
    "--json_body_file",
    required=False,
    help="File containing JSON to use in the body of a POST call.",
)
parser.add_argument(
    "--key_name",
    required=True,
    help="Name of API key - use a different name for each site.",
)
parser.add_argument(
    "--method",
    required=False,
    choices=["GET", "POST"],
    default="GET",
    help="HTTP method",
)
parser.add_argument(
    "--output_file",
    required=False,
    help="Path and filename of saved file. Output sent to STDOUT if not specified.",
)
parser.add_argument(
    "--setup", action="store_true", help="Initial setup to obtain access token."
)
parser.add_argument("--site", required=False, choices=["PubMLST", "Pasteur"])
parser.add_argument(
    "--token_dir",
    required=False,
    default="./.bigsdb_tokens",
    help="Directory into which keys and tokens will be saved.",
)
parser.add_argument("--url", required=False, help="URL for API call.")
args = parser.parse_args()


def main():
    check_required_args(args)
    check_dir(args.token_dir)
    if args.setup:
        (access_token, access_secret) = get_new_access_token()
        if not access_token or not access_secret:
            raise PermissionError("Cannot get new access token.")
    (token, secret) = retrieve_token("session")
    if not token or not secret:
        (token, secret) = get_new_session_token()
    if not args.setup:
        get_route(args.url, token, secret)


def check_required_args(args):
    if args.setup:
        if not args.site:
            parser.error("--site is required for setup")
        if not args.db:
            parser.error("--db is required for setup")
    else:
        if not args.url:
            parser.error("--url is required")
    if args.json_body:
        if args.method != "POST":
            parser.error("You cannot use --json_body with --method=GET")
        if args.json_body_file:
            parser.error("You cannot use both --json_body and --json_body_file")

    if args.json_body_file:
        if args.method != "POST":
            parser.error("You cannot use --json_body_file with --method=GET")


def is_valid_json(json_string):
    try:
        json.loads(json_string)
        return True
    except ValueError:
        return False


def trim_url_args(url):
    if not "?" in url:
        return url, {}
    trimmed_url, param_string = url.split("?")
    params = parse_qs(param_string)
    params = {k: int(v[0]) for k, v in params.items()}

    return trimmed_url, params


def get_route(url, token, secret):
    (client_key, client_secret) = get_client_credentials()
    session = OAuth1Session(
        client_key, client_secret, access_token=token, access_token_secret=secret
    )
    trimmed_url, request_params = trim_url_args(url)
    if args.method == "GET":
        r = session.get(
            trimmed_url,
            params=request_params,
            headers={"User-Agent": "BIGSdb downloader"},
        )
    else:
        json_body = "{}"
        if args.json_body:
            json_body = args.json_body
        elif args.json_body_file:
            with open(args.json_body_file, "r") as file:
                json_body = file.read()
        if not is_valid_json(json_body):
            parser.error("Body does not contain valid JSON")
        r = session.post(
            trimmed_url,
            params=request_params,
            data=json_body,
            headers={
                "Content-Type": "application/json",
                "User-Agent": "BIGSdb downloader",
            },
            header_auth=True,
        )

    if r.status_code == 200 or r.status_code == 201:
        if args.output_file:
            try:
                with open(args.output_file, "w") as file:
                    if re.search("json", r.headers["content-type"], flags=0):
                        file.write(str(r.json()))
                    else:
                        file.write(r.text)
            except IOError as e:
                sys.stderr.write(f"An error occurred while writing to the file: {e}\n")
        else:
            if re.search("json", r.headers["content-type"], flags=0):
                print(r.json())
            else:
                print(r.text)
    elif r.status_code == 400:
        sys.stderr.write("Bad request")
        sys.stderr.write(r.json()["message"])
        sys.exit(1)
    elif r.status_code == 401:
        if re.search("unauthorized", r.json()["message"]):
            sys.stderr.write("Access denied - client is unauthorized\n")
            sys.exit(1)
        else:
            sys.stderr.write(r.json()["message"] + "\n")
            sys.stderr.write("Invalid session token, requesting new one...\n")
            (token, secret) = get_new_session_token()
            get_route(url, token, secret)
    else:
        sys.stderr.write(f"Error: {r.text}\n")
        sys.exit(1)


def check_dir(directory):
    if os.path.isdir(directory):
        if os.access(directory, os.W_OK):
            return
        else:
            raise PermissionError(
                f"The token directory '{directory}' exists but is not writable."
            )
    else:
        try:
            os.makedirs(directory)
            os.chmod(directory, stat.S_IRWXU)  # Set permissions to 0700
        except OSError as e:
            raise PermissionError(
                f"Failed to create token directory '{directory}': {e}"
            )


def retrieve_token(token_type):
    file_path = Path(f"{args.token_dir}/{token_type}_tokens")
    if file_path.is_file():
        config = configparser.ConfigParser(interpolation=None)
        config.read(file_path)
        if config.has_section(args.key_name):
            token = config[args.key_name]["token"]
            secret = config[args.key_name]["secret"]
            return (token, secret)
    return (None, None)


def get_new_session_token():
    file_path = Path(f"{args.token_dir}/session_tokens")
    (access_token, access_secret) = retrieve_token("access")
    if not access_token or not access_secret:
        (access_token, access_secret) = get_new_access_token()
    service = get_service()
    (client_key, client_secret) = get_client_credentials()
    db = get_db_value()
    url = f"{BASE_API[args.site]}/db/{db}/oauth/get_session_token"
    session_request = OAuth1Session(
        client_key,
        client_secret,
        access_token=access_token,
        access_token_secret=access_secret,
    )
    r = session_request.get(url, headers={"User-Agent": "BIGSdb downloader"})
    if r.status_code == 200:
        token = r.json()["oauth_token"]
        secret = r.json()["oauth_token_secret"]
        config = configparser.ConfigParser(interpolation=None)
        if file_path.is_file():
            config.read(file_path)
        config[args.key_name] = {"token": token, "secret": secret}
        with open(file_path, "w") as configfile:
            config.write(configfile)

        return (token, secret)
    else:
        sys.stderr.write(
            "Failed to get new session token. " + r.json()["message"] + "\n"
        )
        if args.cron:
            sys.stderr.write("Run interactively to fix.\n")
        if re.search("verification", r.json()["message"]) or re.search(
            "Invalid access token", r.json()["message"]
        ):
            sys.stderr.write("New access token required - removing old one.\n")
            config = configparser.ConfigParser(interpolation=None)
            file_path = Path(f"{args.token_dir}/access_tokens")
            if file_path.is_file():
                config.read(file_path)
                config.remove_section(args.key_name)
                with open(file_path, "w") as configfile:
                    config.write(configfile)
        sys.exit(1)


def get_service():
    db = get_db_value()
    (client_key, client_secret) = get_client_credentials()
    request_token_url = f"{BASE_API[args.site]}/db/{db}/oauth/get_request_token"
    access_token_url = f"{BASE_API[args.site]}/db/{db}/oauth/get_access_token"
    return OAuth1Service(
        name="BIGSdb_downloader",
        consumer_key=client_key,
        consumer_secret=client_secret,
        request_token_url=request_token_url,
        access_token_url=access_token_url,
        base_url=BASE_API[args.site],
    )


def get_new_access_token():
    if args.cron:
        sys.stderr.write(f"No access token saved for {args.key_name}.\n")
        sys.stderr.write("Run interactively to set.\n")
        sys.exit(1)
    file_path = Path(f"{args.token_dir}/access_tokens")
    (request_token, request_secret) = get_new_request_token()
    db = get_db_value()
    print(
        "Please log in using your user account at "
        f"{BASE_WEB[args.site]}?db={db}&page=authorizeClient&oauth_token={request_token} "
        "using a web browser to obtain a verification code."
    )
    verifier = input("Please enter verification code: ")
    service = get_service()
    r = service.get_raw_access_token(
        request_token,
        request_secret,
        params={"oauth_verifier": verifier},
        headers={"User-Agent": "BIGSdb downloader"},
    )
    if r.status_code == 200:
        token = r.json()["oauth_token"]
        secret = r.json()["oauth_token_secret"]
        file_path = Path(f"{args.token_dir}/access_tokens")
        print("Access Token:        " + token)
        print("Access Token Secret: " + secret + "\n")
        print(
            "This access token will not expire but may be revoked by the \n"
            f"user or the service provider. It will be saved to \n{file_path}."
        )
        config = configparser.ConfigParser(interpolation=None)
        if file_path.is_file():
            config.read(file_path)
        config[args.key_name] = {"token": token, "secret": secret}
        with open(file_path, "w") as configfile:
            config.write(configfile)
        return (token, secret)
    else:
        sys.stderr.write("Failed to get new access token." + r.json()["message"])
        sys.exit(1)


def get_db_value():
    db = None
    if args.db:
        db = args.db
    elif args.url:
        match = re.search(r"/db/([^/]+)", args.url)
        if match:
            db = match.group(1)
        else:
            raise ValueError("No db value found in the URL.")
    return db


def get_new_request_token():
    (client_key, client_secret) = get_client_credentials()
    db = get_db_value()
    service = get_service()

    r = service.get_raw_request_token(
        params={"oauth_callback": "oob"}, headers={"User-Agent": "BIGSdb downloader"}
    )
    if r.status_code == 200:
        token = r.json()["oauth_token"]
        secret = r.json()["oauth_token_secret"]
        return (token, secret)
    else:
        sys.stderr.write("Failed to get new request token." + r.json()["message"])
        sys.exit(1)


def get_client_credentials():
    config = configparser.ConfigParser(interpolation=None)
    file_path = Path(f"{args.token_dir}/client_credentials")
    client_id = None
    client_secret = None
    if file_path.is_file():
        config.read(file_path)
        if config.has_section(args.key_name):
            client_id = config[args.key_name]["client_id"]
            client_secret = config[args.key_name]["client_secret"]
    if not client_id:
        if args.cron:
            sys.stderr.write(f"No client credentials saved for {args.key_name}.\n")
            sys.stderr.write("Run interactively to set.\n")
            sys.exit(1)
        client_id = input("Enter client id: ").strip()
        while len(client_id) != 24:
            print("Client ids are exactly 24 characters long.")
            client_id = input("Enter client id: ").strip()
        client_secret = input("Enter client secret: ").strip()
        while len(client_secret) != 42:
            print("Client secrets are exactly 42 characters long.")
            client_secret = input("Enter client secret: ").strip()

        config[args.key_name] = {"client_id": client_id, "client_secret": client_secret}
        with open(file_path, "w") as configfile:
            config.write(configfile)
    return client_id, client_secret


if __name__ == "__main__":
    main()
