#!/usr/bin/python

def main():

	vx = [[0 for i in range(500)] for j in range(500)]
	vy = [[0 for i in range(500)] for j in range(500)]
	(vx[0])[2] = 5
	print (vx[0])[2]


if __name__ == '__main__':
	main()
