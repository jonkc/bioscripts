import sys
import os

#Version 04... rewrote the algorithm in the comparison part... much faster now. Also adding option to output the total only.
#Version 05... major rewrite so that the script processes data in a stream to minimize ram usage

#cwd = os.getcwd()

def file_path(argv):
	cwd = os.getcwd()
	if os.path.isabs(argv):
		return argv
	else:
		return cwd + "/" + argv

#Creates array of the consensus sites to count cleavage sites
def init_array(fname, dist_5, dist_3):
	arrays = []
	with open(fname) as f:
		next(f)
		array = ['chrID', 'start', 'stop', 'center', 'strand', 'p-value']
		columns = next(f).split('\t')
		odd_len = 0
		if (int(columns[4])-int(columns[3]) + 1) % 2 > 0:
			odd_len = True
			region_len = dist_5 + dist_3 + 1
			for i in range(region_len):
				array.append(-dist_5 + i)
			array.append(region_len)
		elif (int(columns[4])-int(columns[3]) + 1) % 2 == 0:
			odd_len = False
			region_len = dist_5 + dist_3
			for i in range(dist_5):
				array.append(-dist_5 + i)
			for i in range(dist_3):
				array.append(i + 1)
			array.append(region_len)
		arrays.append(array)
		array = []

		start = int(columns[3])
		end = int(columns[4])
		for i in [2, 3, 4, 5, 7]:
			#modifies position values of array (duplicate #1)
			if i == 3:
				if odd_len is True:
					if columns[5] == '+':
						array.append((end - start)//2 + start - dist_5)
						#array.append((end - start) + start)
					elif columns[5] == '-':
						array.append((end - start)//2 + start - dist_3)
					else:
						print "WTF!!! shitty data"
				elif odd_len is False:
					if columns[5] == '+':
						array.append((end - start)//2 + start - dist_5)
					elif columns[5] == '-':
						array.append((end - start)//2 + start - dist_3)
					else:
						print "WTF!!! shitty data #2"
			elif i == 4:
				if odd_len is True:
					if columns[5] == '+':
						array.append((end - start)//2 + start + dist_3)
						#array.append((end - start) + start)
					elif columns[5] == '-':
						array.append((end - start)//2 + start + dist_5)
					else:
						print "WTF!!! shitty data #3"
				elif odd_len is False:
					if columns[5] == '+':
						array.append((end - start)//2 + start + dist_3 - 1)
					elif columns[5] == '-':
						array.append((end - start)//2 + start + dist_5 - 1)
					else:
						print "WTF!!! shitty data #4"

				#Get center position of region of interest
				if odd_len is True:
					if columns[5] == '+':
						array.append(array[2] - dist_3)
					elif columns[5] == '-':
						array.append(array[2] - dist_5)
				elif odd_len is False:
					if columns[5] == '+':
						array.append(((array[2] - dist_3) + (array[1] + dist_5))/2.0)
					if columns[5] == '-':
						array.append(((array[2] - dist_5) + (array[1] + dist_3))/2.0)
			#end modification
			else:
				array.append(columns[i])
		for i in range(region_len):
			array.append(0)
		arrays.append(array)

		for line in f:
			array = []
			if '\t' in line:
				columns = line.split('\t')
				start = int(columns[3])
				end = int(columns[4])
				for i in [2, 3, 4, 5, 7]:
					#modifies position values of array (duplicate #2)
					if i == 3:
						if odd_len is True:
							if columns[5] == '+':
								array.append((end - start)//2 + start - dist_5)
								#array.append((end - start) + start)
							elif columns[5] == '-':
								array.append((end - start)//2 + start - dist_3)
							else:
								print "WTF!!! shitty data"
						elif odd_len is False:
							if columns[5] == '+':
								array.append((end - start)//2 + start - dist_5)
							elif columns[5] == '-':
								array.append((end - start)//2 + start - dist_3)
							else:
								print "WTF!!! shitty data #2"
					elif i == 4:
						if odd_len is True:
							if columns[5] == '+':
								array.append((end - start)//2 + start + dist_3)
								#array.append((end - start) + start)
							elif columns[5] == '-':
								array.append((end - start)//2 + start + dist_5)
							else:
								print "WTF!!! shitty data #3"
						elif odd_len is False:
							if columns[5] == '+':
								array.append((end - start)//2 + start + dist_3 - 1)
							elif columns[5] == '-':
								array.append((end - start)//2 + start + dist_5 - 1)
							else:
								print "WTF!!! shitty data #4"

						#Get center position of region of interest
						if odd_len is True:
							if columns[5] == '+':
								array.append(array[2] - dist_3)
							elif columns[5] == '-':
								array.append(array[2] - dist_5)
						elif odd_len is False:
							if columns[5] == '+':
								array.append(((array[2] - dist_3) + (array[1] + dist_5))/2.0)
							if columns[5] == '-':
								array.append(((array[2] - dist_5) + (array[1] + dist_3))/2.0)
					#end modification
					else:
						array.append(columns[i])
				for i in range(region_len):
					array.append(0)
				arrays.append(array)
		return arrays
			#print columns[2]

def count_fragments(infile,fimofile, dist_5, dist_3):
	#array to record counting results
	arrays = init_array(fimofile, dist_5, dist_3)
	region_len = arrays[0][-1]
	regions = len(arrays) - 1

	#get element index for center position
	for i in range(region_len):
		if arrays[0][i + 6] == 0:
			center = i + 6
			break
		elif arrays[0][i + 6] == 1:
			center = i + 6 - 0.5
			break

	#create a list for total, which will be appended to the *arrays* list later
	total = ['SUM', 'N/A', 'N/A','N/A','N/A','N/A']
	for i in range(region_len):
		total.append(0)	

	#comparison part
	current_index = 1
	with open(infile) as f:
		next(f)
		#lines = f.readlines()
		for line in f:
			columns = line.split('\t')

			#if the chromosome don't match, move the current_index for fimo arrays up
			while True:
				if columns[0] != arrays[current_index][0]:
					current_index += 1
					if current_index > regions:
						arrays.append(total)
						return arrays
				else:
					break

			#if the f query position is larger than the current position range of the fimo entry, move index up
			while True:
				if int(columns[1]) > arrays[current_index][2]:
					current_index += 1
					if current_index > regions:
						arrays.append(total)
						return arrays
				else:
					break

			#if in range of fimo entry
			if int(columns[1]) >= arrays[current_index][1] and int(columns[1]) <= arrays[current_index][2]:
				distance_from_center = arrays[current_index][3] - int(columns[1])
				if (abs(distance_from_center) % 1) == 0:
					if arrays[current_index][4] == '+':
						arrays[current_index][center - distance_from_center] += abs(int(columns[2]))
						total[center - distance_from_center] += abs(int(columns[2]))
					elif arrays[current_index][4] == '-':
						arrays[current_index][center + distance_from_center] += abs(int(columns[2]))
						total[center + distance_from_center] += abs(int(columns[2]))
				if (abs(distance_from_center) % 1) > 0:
					if arrays[current_index][4] == '+':
						arrays[current_index][int(center - distance_from_center)] += abs(int(columns[2]))
						total[int(center - distance_from_center)] += abs(int(columns[2]))
					elif arrays[current_index][4] == '-':
						arrays[current_index][int(center + distance_from_center)] += abs(int(columns[2]))
						total[int(center + distance_from_center)] += abs(int(columns[2]))


			''' Inefficient version.
			for i in range(regions):
				if columns[0] == arrays[i + 1][0]:
					if int(columns[1]) > arrays[i + 1][1] and int(columns[1]) < arrays[i + 1][2]:
						#if consensus was found in reverse strand, put the cleavage site in the right array element
						distance_from_center = arrays[i + 1][3] - int(columns[1])
						if (abs(distance_from_center) % 1) == 0:
							if arrays[i + 1][4] == '+':
								arrays[i + 1][center - distance_from_center] += abs(int(columns[2]))
							elif arrays[i + 1][4] == '-':
								arrays[i + 1][center + distance_from_center] += abs(int(columns[2]))
						if (abs(distance_from_center) % 1) > 0:
							if arrays[i + 1][4] == '+':
								arrays[i + 1][center - distance_from_center] += abs(int(columns[2]))
							elif arrays[i + 1][4] == '+':
								arrays[i + 1][center + distance_from_center] += abs(int(columns[2]))

	#Get sum, obsolete. incorporated to the comparison part
	array = ['SUM', 'N/A', 'N/A','N/A','N/A','N/A']
	for x in range(region_len):
		z = 0
		for y in range(regions):
			z += arrays[y + 1][x + 6]
		array.append(z)
	arrays.append(array)
	'''
	arrays.append(total)
	return arrays

def fimo_entry(line, dist_5, dist_3, odd_len, region_len):
	#processes an entry (line) from the fimo output file and returns it as a list
	if line == -1:
		return -1
	array = []
	if '\t' in line:
		columns = line.split('\t')
		start = int(columns[3])
		end = int(columns[4])
		for i in [2, 3, 4, 5, 7]:
			#modifies position values of array (duplicate #2)
			if i == 3:
				if odd_len is True:
					if columns[5] == '+':
						array.append((end - start)//2 + start - dist_5)
						#array.append((end - start) + start)
					elif columns[5] == '-':
						array.append((end - start)//2 + start - dist_3)
					else:
						print "WTF!!! shitty data"
				elif odd_len is False:
					if columns[5] == '+':
						array.append((end - start)//2 + start - dist_5)
					elif columns[5] == '-':
						array.append((end - start)//2 + start - dist_3)
					else:
						print "WTF!!! shitty data #2"
			elif i == 4:
				if odd_len is True:
					if columns[5] == '+':
						array.append((end - start)//2 + start + dist_3)
						#array.append((end - start) + start)
					elif columns[5] == '-':
						array.append((end - start)//2 + start + dist_5)
					else:
						print "WTF!!! shitty data #3"
				elif odd_len is False:
					if columns[5] == '+':
						array.append((end - start)//2 + start + dist_3 - 1)
					elif columns[5] == '-':
						array.append((end - start)//2 + start + dist_5 - 1)
					else:
						print "WTF!!! shitty data #4"

				#Get center position of region of interest
				if odd_len is True:
					if columns[5] == '+':
						array.append(array[2] - dist_3)
					elif columns[5] == '-':
						array.append(array[2] - dist_5)
				elif odd_len is False:
					if columns[5] == '+':
						array.append(((array[2] - dist_3) + (array[1] + dist_5))/2.0)
					if columns[5] == '-':
						array.append(((array[2] - dist_5) + (array[1] + dist_3))/2.0)
			#end modification
			else:
				array.append(columns[i])
		for i in range(region_len):
			array.append(0)
	return array

def stream_version(iname, fname, dist_5, dist_3):
	header=[]
	temp=[]
	temps=[]

	#create header
	with open(fname) as f:
		next(f)
		header = ['chrID', 'start', 'stop', 'center', 'strand', 'p-value']
		columns = next(f).split('\t')
		odd_len = 0
		if (int(columns[4])-int(columns[3]) + 1) % 2 > 0:
			odd_len = True
			region_len = dist_5 + dist_3 + 1
			for i in range(region_len):
				header.append(-dist_5 + i)
			header.append(region_len)
		elif (int(columns[4])-int(columns[3]) + 1) % 2 == 0:
			odd_len = False
			region_len = dist_5 + dist_3
			for i in range(dist_5):
				header.append(-dist_5 + i)
			for i in range(dist_3):
				header.append(i + 1)
			header.append(region_len)
	#get index for center position
	for i in range(region_len):
		if header[i + 6] == 0:
			center = i + 6
			break
		elif header[i + 6] == 1:
			center = i + 6 - 0.5
			break
	print header

	with open(iname) as z, open(fname) as f:
		next(z)
		next(f)
		fimo_list = []
		#get first few fimo entries
		fimo_list.append(fimo_entry(next(f), dist_5, dist_3, odd_len, region_len))
		while True:
			fimo_entry = fimo_entry(next(f), dist_5, dist_3, odd_len, region_len)
			if fimo_entry != -1:
				fimo_list.append(fimo_entry)
			else:
				break
			if fimo_list[(len(fimo_list)-1)][0] != fimo_list[0][0]:
				break
			elif fimo_list[(len(fimo_list)-1)][1] > fimo_list[0][2]:
				break

		for line in z:
			query = line.split('\t')
			if int(query[1]) > fimo_list[0]:
				print fimo_list[0]
				fimo_list = fimo_list[1:]

				while True:
					fimo_entry = fimo_entry(next(f), dist_5, dist_3, odd_len, region_len)
					if fimo_entry != -1:
						fimo_list.append(fimo_entry)
					else:
						break
					if fimo_list[(len(fimo_list)-1)][0] != fimo_list[0][0]:
						break
					elif fimo_list[(len(fimo_list)-1)][1] > fimo_list[0][2]:
						break


			for entry in fimo_list:
				#if chrid dont match move on to next line in z
				if query[0] != entry[0]:
					break
				#if query pos is < than 5post, move on
				if int(query[1]) < entry[1]:
					break


		'''
		while True:


		#current_query = fimo_query(next(f), dist_5, dist_3, odd_len, region_len)
		next_entry = fimo_entry(next(f), dist_5, dist_3, odd_len, region_len)
		while True:
			current_entry = next_entry
			next_entry = fimo_entry(next(f, -1), dist_5, dist_3, odd_len, region_len)
			#print current_query

			#cycle through queries list
			if len(queries) > 0:
				for i in range(len(queries)):
					if queries[i][0] == current_entry[0]:

			while True:
				if len(queries) > 0:
					for q in queries:
						if q[0] == current_query:

						else
				query = next(z).split('\t')
				if query[0] != 


			if next_query == -1:
				break
			'''


in_file = file_path(sys.argv[1])
fimo_file = file_path(sys.argv[2])
len_5 = int(sys.argv[3])
len_3 = int(sys.argv[4])
if len(sys.argv) > 5:
	if sys.argv[5] == 'sum-only':
		sum_only = True
	else:
		sum_only = False
else:
	sum_only = False

stream_version(in_file, fimo_file, len_5, len_3)

'''
arrays = count_fragments(in_file, fimo_file, len_5, len_3)

#print results
if sum_only is True:
	print('\t'.join(map(str,arrays[0])))
	print('\t'.join(map(str,arrays[len(arrays)-1])))
else:
	for array in arrays:
		print('\t'.join(map(str,array)))
'''