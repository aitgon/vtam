def insert_table(session, model, obj):
	try:  # checks if exists Phenotype in db
		session.query(model).filter_by(**obj).one()
	except:  # if not, add
		session.add(model(**obj))


def fasta_writer(session, model, file_name):
	name = file_name.replace(".fasta", "_count.fasta")
	file = open(name, "w")
	j = 1
	line = []
	for element in session.query(model).all():
		file.write("> Read number " + str(j) + " count: " + str(element.count))
		file.write("\n")
		# for i in range(len(element.sequence)):
		# 	if (i%80) == 0:
		# 		file.write(sequence)
		# 		file.write('\n')
		file.write(element.sequence + '\n')
		j += 1