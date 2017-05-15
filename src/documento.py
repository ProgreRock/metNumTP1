# SIMPLE
# Plots square root function from 0 to 100
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pylab
import sys


def simple():

	if len(sys.argv) <2:
		sys.stderr.write('hay que pasar un archivo como parametro\n')
		sys.stderr.write('el formato del archivo debe ser:\nx_0 y_0\nx_1 y_1\n...\n')
		
		sys.exit(1)
	
	else:
		with open(sys.argv[1],'r') as file:
			data = file.readlines()
	
		for line in data:
			palabra = line.split()
			matplotlib.pyplot.scatter(palabra[0], palabra[1])
	
		

		axes = matplotlib.pyplot.gca()
		axes.set_xlabel('Cantidad de corridas')
		axes.set_ylabel('Tiempo Corridas')
		pylab.savefig("Mediciones.png")
	

simple()

