import sqlite3

__all__ = ["StationDB",]

def _single_out(array):
    """ If an array has one element, pull it out.

        Parameters
        ----------
        array: array-like
            Input array

        Returns
        -------
        element or array
            None if the array has no element.
            When the array has one element, the element is returned.
            When the array has more then one element, just return
            the array back.
    """
    if len(array) < 1:
        return None
    if len(array) == 1:
        return array[0]
    else:
        return array

class StationDB(object):
    """ Read a station DB in sqlite format, and provide station information
    """
    def __init__(self, fname):
        """ Constructor

            fname:
                sqlite3 station db file
        """
        # uri = "file:%s?mode=ro" % fname
        # self._conn = sqlite3.connect(uri, uri=True)
        self._conn = sqlite3.connect(fname)

    def get_alias(self, station_id):
        """ Get an alias of a station with an ID
            station_id:
                station name
        """
        c = self._conn.cursor()
        sql = "SELECT alias_id FROM stations_utm WHERE Id ='%s'" % station_id
        c.execute(sql)
        results = [row[0] for row in c.fetchall()]
        # return results
        return _single_out(results)

    def get_long_name(self, station_id):
        """ Get a long name of a station with an ID

            Parameters
            ----------
            station_id:
                station name

            Returns
            -------
            str
                a long name of a station
        """
        c = self._conn.cursor()
        sql = "SELECT Name FROM stations_utm WHERE Id ='%s'" % station_id
        c.execute(sql)
        results = [row[0] for row in c.fetchall()]
        return _single_out(results)

    def get_station_id_from_alias(self, alias):
        """
        """
        c = self._conn.cursor()
        sql = "SELECT Id FROM stations_utm WHERE alias_id ='%s'" % alias
        c.execute(sql)
        results = [row[0] for row in c.fetchall()]
        # return _single_out(results)
        return results

    def __del__(self):
        self._conn.close()
