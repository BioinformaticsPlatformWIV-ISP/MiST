from mist.app.dbs.bigsdbauthdownloader import BIGSDbAuthDownloader
from mist.app.dbs.bigsdbdownloader import BIGSDbDownloader
from mist.app.dbs.cgmlstorgdownloader import CgMLSTOrgDownloader
from mist.app.dbs.enterobasedownloader import EnteroBaseDownloader

DOWNLOADERS = {
    'cgmlstorg': CgMLSTOrgDownloader,
    'enterobase': EnteroBaseDownloader,
    'bigsdb': BIGSDbDownloader,
    'bigsdb_auth': BIGSDbAuthDownloader,
}

BASE_URL_BY_DB = {
    'cgmlstorg': {'base': 'https://www.cgmlst.org/ncs'},
    'enterobase': {'base': 'https://enterobase.warwick.ac.uk/schemes'},
    'bigsdb': {'pubmlst': 'https://rest.pubmlst.org', 'pasteur': 'https://bigsdb.pasteur.fr/api'},
    'bigsdb_auth': {'pubmlst': 'https://rest.pubmlst.org', 'pasteur': 'https://bigsdb.pasteur.fr/api'},
}
