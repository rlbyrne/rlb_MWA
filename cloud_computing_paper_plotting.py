#!/usr/bin/python

import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt


def plot_iops():

    path = '/Users/ruby/Downloads/CloudWatch_Data'

    read_ops_file = 'volume_read_ops_4pol_decon.json'
    write_ops_file = 'volume_write_ops_4pol_decon.json'

    read_data = pd.read_json(f'{path}/{read_ops_file}')
    read_timestamps = np.array([
        read_data['Datapoints'][ind]['Timestamp']
        for ind in range(len(read_data['Datapoints']))
    ])
    read_values = np.array([
        read_data['Datapoints'][ind]['Sum']
        for ind in range(len(read_data['Datapoints']))
    ])
    print(np.shape(read_values))
    read_times = np.array([datetime.strptime(
        read_timestamps[ind], '%Y-%m-%dT%H:%M:%SZ'
    ) for ind in range(len(read_timestamps))])
    read_time_diff = read_times - min(read_times)
    sort_inds = np.argsort(read_time_diff)
    times_sorted = np.array([
        read_time_diff[ind].total_seconds() for ind in sort_inds
    ])
    read_values = read_values[sort_inds]

    write_data = pd.read_json(f'{path}/{read_ops_file}')
    write_timestamps = np.array([
        read_data['Datapoints'][ind]['Timestamp']
        for ind in range(len(read_data['Datapoints']))
    ])
    write_values = np.array([
        read_data['Datapoints'][ind]['Sum']
        for ind in range(len(read_data['Datapoints']))
    ])
    write_times = np.array([datetime.strptime(
        write_timestamps[ind], '%Y-%m-%dT%H:%M:%SZ'
    ) for ind in range(len(read_timestamps))])
    write_time_diff = write_times - min(write_times)
    sort_inds = np.argsort(write_time_diff)
    write_values = write_values[sort_inds]

    total_vals = read_values+write_values
    time_interval = float(max(times_sorted))/len(times_sorted)
    iops = total_vals/time_interval
    iops = np.append(np.array([0]), iops)
    times_sorted = np.append(times_sorted, np.array([times_sorted[-1]+time_interval]))
    print(time_interval)

    plt.figure()
    plt.plot(times_sorted/60., iops, "-o", markersize=1)
    plt.xlabel("time (m)")
    plt.ylabel("IOPS")
    plt.show()


def plot_cpu_usage():

    path = '/Users/ruby/Downloads/CloudWatch_Data'
    cpu_utilization_file = 'cpu_utilization_4pol_decon.json'

    data = pd.read_json(f'{path}/{cpu_utilization_file}')
    timestamps = np.array([
        data['Datapoints'][ind]['Timestamp']
        for ind in range(len(data['Datapoints']))
    ])
    values = np.array([
        data['Datapoints'][ind]['Maximum']
        for ind in range(len(data['Datapoints']))
    ])
    print(np.shape(values))
    times = np.array([datetime.strptime(
        timestamps[ind], '%Y-%m-%dT%H:%M:%SZ'
    ) for ind in range(len(timestamps))])
    time_diff = times - min(times)
    sort_inds = np.argsort(time_diff)
    times_sorted = np.array([
        time_diff[ind].total_seconds() for ind in sort_inds
    ])
    values = values[sort_inds]

    plt.figure()
    plt.plot(values, "-o", markersize=1)
    plt.xlabel("time (m)")
    plt.ylabel("CPU usage (%)")
    plt.show()


def plot_ram_usage():

    path = '/Users/ruby/Downloads/CloudWatch_Data'
    ram_use_file = '1131731632_ram_usage_365_3.84.16.70.txt'

    file = open(f'{path}/{ram_use_file}', 'r')
    file_data = file.readlines()
    file.close()

    'Fri Feb  7 09:14:24 UTC 2020'
    times = np.array([line.rstrip() for line in file_data if 'UTC' in line])
    times = np.array([datetime.strptime(
        times[ind], '%a %b %d %H:%M:%S UTC %Y'
    ) for ind in range(len(times))])
    data = np.array([
        float(line.rstrip()) for line in file_data
        if 'UTC' not in line and 'RAM' not in line
    ])

    time_diff = times - min(times)
    sort_inds = np.argsort(time_diff)
    times_sorted = np.array([
        time_diff[ind].total_seconds() for ind in sort_inds
    ])
    data_sorted = data[sort_inds]
    plt.figure()
    plt.plot(times_sorted/60., data_sorted/1024., "-o", markersize=1)
    plt.xlabel("time (m)")
    plt.ylabel("RAM usage (GB)")
    plt.show()
    plt.close()


if __name__ == '__main__':
    plot_ram_usage()
    