# Licensed under GPL v3. See license file for details.
'''
This module contains classes that wrap astroquery queries.

For each catalog or service (e.g. Simbad, the 4XMM catalog, ...) we define one class.
Each of those derives from a base class called `AstroService` with three main functionalities:

- query on web or, if present, load result of existing query from disk
- produce short, meaningful description of query result
- add symbols to a plot

For each query, we derive from that base class and change a few details, e.g. the plotting symbol.
Simple changes can be made by a simple declaration by just defining a few class attributes, more
complex changes can be done by overriding the query or plot method.
'''
import os
import functools
import datetime

from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord

from astroquery.simbad import Simbad as aSimbad
from astroquery.gemini import Observations as Gemini
from astroquery.alma import Alma
from astroquery.nrao import Nrao
from astroquery.mast import Observations as Mast
from astroquery.casda import Casda  #  Australian Square Kilometre Array Pathfinder (ASKAP)
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT == -1
from astroquery.utils.commons import TableList
from astroquery.eso import Eso

markerdefaults = {'s': 100, 'facecolor': 'none'}


class TableNotLoadedError(Exception):
    '''No data table has been obtained through query or loading yet'''
    pass


def table_present(func):
    '''Decorator to check if a table has been loaded'''
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self.table is None:
            raise TableNotLoadedError(f'Data table for {self.name} not yet loaded.')
        return func(self, *args, **kwargs)
    return wrapper


def skip_if_empty(func):
    '''Decorator to skip a function if the table of sources is empty'''
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if len(self.table) == 0:
            print(f"Zero lines in {self.name} data table!")
            return
        else:
            return func(self, *args, **kwargs)
    return wrapper


class AstroService():
    name = None
    'Human readable name for this dataset. If `None` use the class name.'

    radius = 2 * u.arcmin
    'Query radius. The default is set to be a litte bigger than the FOV.'

    markerargs = {}
    'Keyword arguments passed to plt.scatter to plot symbols'

    table = None
    'Data table'

    @property
    @table_present
    def ra(self):
        'Get names of first column that looks like a RA'
        for col in self.table.columns:
            if 'ra' in col or 'RA' in col:
                return col
        else:
            raise KeyError('No column looks like a RA.')

    @property
    @table_present
    def dec(self):
        'Get name of first column that looks like a DEC'
        for col in self.table.columns:
            if 'de' in col or 'DE' in col:
                return col
        else:
            raise KeyError('No column looks like a DE.')

    def __repr__(self):
        if self.table is None:
            return f'{self.__class__}: No table loaded'
        else:
            return self.table.__repr__()

    @table_present
    def description(self):
        'Return string with a very short description of the data table'
        return f' {len(self.table)}\n'

    def __init__(self, name=None):
#        print(self.name)
        if self.name is None:
            self.name = self.service.__class__.__name__
            
    @table_present
    def write(self, directory):
        'Write data table to directory. File name and format are fixed to name.ecsv'
        self.table.meta['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d")
        self.table.write(os.path.join(directory, self.name + '.ecsv'))

    def read(self, directory):
        'Read data table from directory'
        # In case any error happens when reading, we want to make sure it's 
        # at least set to None and we don't accidentially continue with an old table
        self.table = None
        self.table = Table.read(os.path.join(directory, self.name + '.ecsv'))

    def query(self, center):
        'Run astroquery around center with the radius set in self.radius'
        try:
            # Most queries return a TableList (e.g. Vizier query)
            # but some return a single table (e.g. Simbad query)
            out = self.service.query_region(center, radius=self.radius)
            #if type(self.table) == TableList and len(self.table) == 1:
            #    self.table = self.table[0]

            if out is None:
                # Astroqueries way of saying "no match"
                self.table = Table()
            else:
                self.table = out
        except Exception as e:
            print(e)
            # Query failed. What to do?
            self.table = Table()

    def read_or_query(self, directory, center):
        '''Read file with data is present or perform astroquery is not'''
        try:
            self.read(directory)
        except (FileNotFoundError, StopIteration) as e:
            if type(e) == FileNotFoundError:
                self.query(center)
                self.write(directory)
            else:
                self.table = Table()

    @table_present
    @skip_if_empty
    def plot(self, ax, **kwargs):
        kwargs = markerdefaults | self.markerargs | kwargs
        ax.scatter(self.table[self.ra], self.table[self.dec], 
                   transform=ax.get_transform('icrs'), label=self.name, **kwargs)

class VCatalog(AstroService):
    service = Vizier
    def query(self, center):
        self.table = Vizier.query_region(center, width=self.radius, catalog=self.catalog, ROW_LIMIT=-1)
        if len(self.table) == 0:
            # Empty table List -- no match
            self.table = Table()
        elif len(self.table) == 1:
            self.table = self.table[0]
        else:
            raise Exception("No idea what this means... Check query!")


class MAST(AstroService):

    service = Mast
    name = 'MAST' # already initiated in the __init__ function in astroservice
    markerargs = {'s': 100, 'edgecolor': (.5, 0, 0), 'marker': 'o'}

    @table_present
    @skip_if_empty
    def description(self):

        table = self.table  # just save some typing in the lines below...
        out = ''
        if 'TESS' in table['obs_collection']:
            out +='  TESS: yes\n'
        else:
            out +='  TESS: no\n'
        
        if 'GALEX' in table['obs_collection']:
            out += '  GALEX: yes\n'
        else:
            out += '  GALEX: no\n'
        
        table = table[table['obs_collection'] != 'TESS']
        table = table[table['obs_collection'] != 'GALEX']
        ### Remove services from our two HST program IDs here!!!
        table = table[~((table['proposal_id'] == '15888') & (table['obs_collection'] == 'HST'))]
        table = table[~((table['proposal_id'] == '16359') & (table['obs_collection'] == 'HST'))]
        table_grouped = table.group_by('obs_collection')
        for key, group in zip(table_grouped.groups.keys, table_grouped.groups):
            out += f'  {key["obs_collection"]}: {len(group)}\n'
        return out

#    @property
#    @table_present
#    def plot(self, ax, **kwargs):
#        kwargs = markerdefaults | self.markerargs | kwargs
#        trans = ax.get_transform('icrs')
#        for row in self.table:
#            print(row['s_ra'], row['s_dec'])
#            print(len(row))
#            print(row.info())
#            coord = SkyCoord(row['s_ra'] + ' ' + row['s_dec'], unit=(u.hourangle, u.deg))
#        ax.scatter(self.table['s_ra'], self.table['s_dec'], c='k',
#                   transform=trans, label=self.name, **kwargs)


class Simbad(AstroService):
    service = aSimbad
    name = 'Simbad'
    ra = 'RA_d'  # not in default Simbad return values, but added with 'coo(d)1 in __init__
    dec = 'DEC_d'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Overwriting class level variable "service" with an actual
        # instance of that class that is specific to this object
        self.service = self.service()
        self.service.add_votable_fields('bibcodelist')
        self.service.add_votable_fields('coo(d)')

    @table_present
    def description(self):
        return f'  objects: {len(self.table)}\n'

    @table_present
    @skip_if_empty
    def plot(self, ax, **kwargs):
        trans = ax.get_transform('icrs')
        for row in self.table:
            coord = SkyCoord(row['RA'] + ' ' + row['DEC'], unit=(u.hourangle, u.deg))
            name = row['MAIN_ID']
            if (row['BIBLIST'] < 50) & (len(name) > 10):
                name = name[:5] + '...'
            ax.text(coord.ra.deg, coord.dec.deg, name, transform=trans, fontsize=12)
        
        # Three catagories with hardcoded boundaries.
        # Could be made more general, but for just three is the easist to just copy and paste
        ind = self.table['BIBLIST'] < 10
        if ind.sum() > 0:
            ax.scatter(self.table['RA_d'][ind], self.table['DEC_d'][ind], c='k', transform=trans,
                    facecolors='none', s=25, 
                    label=f'{self.name} (n < 10)')
        
        ind = (self.table['BIBLIST'] >= 10) & (self.table['BIBLIST'] < 50)
        if ind.sum() > 0:
            ax.scatter(self.table['RA_d'][ind], self.table['DEC_d'][ind], c='k', transform=trans,
                    facecolors='none', s=75, 
                    label=f'{self.name} (10 < n < 50)')

        ind = self.table['BIBLIST'] > 50
        if ind.sum() > 0:
            ax.scatter(self.table['RA_d'][ind], self.table['DEC_d'][ind], c='k', transform=trans,
                    facecolors='none', s=150, 
                    label=f'{self.name} (n > 50)')
                    

class ESO(AstroService):
    service = Eso()
    service.login("abinks", store_password=True)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Overwriting class level variable "service" with an actual
        # instance of that class that is specific to this object
        self.service = self.service()
        self.mission = self.mission
        self.name = self.name

    def query(self, center):
        self.table = self.service.query_instrument(self.mission,coord1=center.ra.deg, coord2=center.dec.deg, box=2./60)
        if self.table is None:
        # Empty table List -- no match
            self.table = Table()
        elif self.table is not None:
            self.table = self.table
        else:
            raise Exception("No idea what this means... Check query!")
            
class RXS2(VCatalog):
    catalog = 'J/A+A/588/A103/cat2rxs'
    name = '2RXS'  # differs from class name, because class names can't start with numbers
    markerargs = {'edgecolor':'xkcd:azure', 'marker': '^'}

class XMM4(VCatalog):
    catalog = 'IX/59'
    name = '4XMM'  # differs from class name, because class names can't start with numbers
    markerargs = {'edgecolor':'xkcd:crimson', 'marker': '^'}

class CSC2(VCatalog):
    catalog = 'IX/57/csc2master'
    name = 'CSC2'
    markerargs = {'edgecolor':'xkcd:dark green', 'marker': '^'}

class IPHAS(VCatalog):
    catalog = 'II/321/iphas2'
    name = 'IPHAS'
    markerargs = {'edgecolor': 'xkcd:azure', 'marker': '*'}
    
class GaiaEDR3(VCatalog):
    catalog = 'I/350/gaiaedr3'
    name = 'GaiaEDR3'
    markerargs = {'edgecolor': 'xkcd:crimson', 'marker': '*'}
    
class giraffe(ESO):
    mission = 'giraffe'
    name = 'ESO_giraffe'
    markerargs = {'edgecolor': 'xkcd:azure', 'marker': 'D'}

class uves(ESO):
    mission = 'uves'
    name = 'ESO_uves'
    markerargs = {'edgecolor': 'xkcd:crimson', 'marker': 'D'}

class xshooter(ESO):
    mission = 'xshooter'
    name = 'ESO_xshooter'
    markerargs = {'edgecolor': 'xkcd:dark green', 'marker': 'D'}
    
class ALMA(AstroService):
    service = Alma
    name = 'ALMA'
    markerargs = {'edgecolor': 'xkcd:azure', 'marker': 'H'}

