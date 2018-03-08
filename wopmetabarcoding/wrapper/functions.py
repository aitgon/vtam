from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


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


def insert_read(session, model, file):
	with open(file, 'r') as file_fasta:
		i = 1
		sequence = ""
		liste_reads = []
		for line in file_fasta:
			if ">" in line:
				if i == 1:
					sequence_id = line.strip()
					sequence_id = sequence_id.replace('>', '')
					continue
				else:
					liste_reads.append({'sequence_id': sequence_id, 'sequence': sequence})
					sequence = ""
					sequence_id = line.strip()
					sequence_id = sequence_id.replace('>', '')
			else:
				seq = line.replace('\n', '')
				sequence += seq
			i += 1
	liste_reads.append({'sequence_id': sequence_id, 'sequence': sequence})
	for dico in liste_reads:
		for element in session.query(model).all():
			if element.sequence_id in dico.get("sequence_id"):
				session.query(model).filter(model.sequence_id == element.sequence_id).update({model.read: dico.get("sequence")})


def trim_reads(session, model):
	for element in session.query(model).all():
		read_sequence = element.read
		target = element.qrow
		trimmed_read = read_sequence.replace(target, '')
		session.query(model).filter(model.sequence_id == element.sequence_id).update({model.read: trimmed_read})


def reverse_complement(session, model):
	for element in session.query(model).all():
		my_read = Seq(element.read, IUPAC.ambiguous_dna)
		reverse_read = my_read.reverse_complement()
		session.query(model).filter(model.sequence_id == element.sequence_id).update({model.reverse: str(reverse_read)})


def create_forward_fasta(session, model, file):
	with open(file, 'w') as fasta_file:
		for line in session.query(model).all():
			fasta_file.write(">" + line.sequence_id + "\n")
			fasta_file.write(line.reverse)
			fasta_file.write("\n")