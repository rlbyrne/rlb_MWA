#!/usr/bin/python

import os


def main():

    obsfile_name = '/nfs/eor-00/h1/rbyrne/sidelobe_survey_obsIDs.txt'
    obsfile = open(obsfile_name, 'r')
    obsids = [line.split(' ') for line in obsfile.readlines()]
    obsids = [obs[0] for obs in obsids]
    obsfile.close()

    # t = Time([int(obsid) for obsid in obsids], format='gps', scale='utc')
    # jdates = t.jd

    for obs in obsids:
        os.system('python get_obs_info.py --start {} --limit 1 -v '
                  '--full'.format(obs))


if __name__ == '__main__':
    main()
