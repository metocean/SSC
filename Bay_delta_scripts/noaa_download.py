# Script to download NOAA water level data

import urllib2
import bs4
import calendar

def retrieve_csv(url):
    response = urllib2.urlopen(url)
    return response.read()

def retrieve_table(url):
    done = False
    while not done:
        try:
            soup = bs4.BeautifulSoup(urllib2.urlopen(url))
            done = True
        except urllib2.URLError:
            print "Failed to retrieve %s" % url
            print "Try again..."
            pass

    table = soup.table.pre.pre.string
    return table

def write_table(table, fname, first):
    f = open(fname, 'a')
    # Remove the Error line
    table = table[:table.find("Error")]
    if table[-1] != '\n':
        table += '\n'
    if first:
        f.write(table)
    else:
        pos = table.find('\n')
        f.write(table[pos+1:])
    f.flush()
    f.close()

def write_header(fname, headers):
    f = open(fname, 'w')
    for key, value in headers.iteritems():
        buf = "# \"%s\"=\"%s\"\n" % (key, value)
        f.write(buf)
    f.flush()
    f.close()

def retrieve_data(years, id_):
    """ Download stage data from NOAA tidesandcurrents, and save it to
        NOAA CSV format.

        Parameters
        ----------
        years: array-like
            years to retrieve
        id_:
            station id
    """
    fname = "%s_gageheight.txt" % id_
    first = True
    headers = {"agency": "noaa",
               "unit": "meter",
               "datum": "NAVD",
               "station_id": "%s" % id_,
               "item": "elev",
               "timezone": "LST",
               "source": "http://tidesandcurrents.noaa.gov/"
              }
    for year in years:
        for month in range(1, 13):
            month_range = calendar.monthrange(year, month)
            date_start = "%4d%02d01" % (year, month)
            date_end = "%4d%02d%02d" % (year, month, month_range[1])

            datum = "NAVD"
            url = "http://tidesandcurrents.noaa.gov/api/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&station=%s&begin_date=%s&end_date=%s&datum=%s&units=metric&time_zone=LST&format=csv" % (id_, date_start, date_end, datum)
            print "Retrieving %s, %s, %s..." % (id_, date_start, date_end)
            print "URL:", url
            # raw_table = retrieve_table(url)
            raw_table = retrieve_csv(url)
            if raw_table[0] == '\n':
                datum = "STND"
                url = "http://tidesandcurrents.noaa.gov/api/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&station=%s&begin_date=%s&end_date=%s&datum=%s&units=metric&time_zone=LST&format=csv" % (id_, date_start, date_end, datum)
                print "Retrieving Station %s, from %s to %s..." % (id_, date_start, date_end)
                print "URL:", url
                # raw_table = retrieve_table(url)
                raw_table = retrieve_csv(url)
            print "Done retrieving."

            if first:
                headers["datum"] = datum
                write_header(fname, headers)
            write_table(raw_table, fname, first)
            first = False

if __name__ == "__main__":
    stage_stations = [
        "9414290",  # San Francisco
        "9414750",  # Alameda
        "9414523",  # Redwood city
        "9414575",  # Coyote Creek
        "9414863",  # Richmond
        "9415144",  # Port Chicago
        "9415102",  # Martinez-Amorco
        ]

    years = (2013, 2014)
    for id_ in stage_stations:
        retrieve_data(years, id_)
