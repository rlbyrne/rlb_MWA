#!/usr/bin/python

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ao_model_reader as ao_read

def main():

	name = 'PicA-88comp'
	source_models = ao_read.read_model(name)
	
	print source_models[0].get('components')[0]
	
if __name__ == '__main__':
	main()
