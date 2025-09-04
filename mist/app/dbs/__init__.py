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
