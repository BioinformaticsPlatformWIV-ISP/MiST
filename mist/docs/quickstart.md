# Quickstart

## 1. Download an MLST scheme

Download from EnteroBase.

```
mkdir yersinia_mlst
cd yersinia_mlst
mmst download --origin enterobase --url {URL} --include-profiles
```

After this you should have a file called `fasta_list.txt` containing the paths to the downloaded FASTA file.

**Note:** Scripts for downloading schemes form other sources are available [here](LINK).

## 2. Create indices

{TOOL_NAME} requires a pre-indexed database. This can be done using the following command:

```
mmst index --fasta-list fasta_list.txt
```

## 3. Query the scheme

A scheme can be queried using the following command:

```

```

## 4. Output

The default output is in JSON format.
