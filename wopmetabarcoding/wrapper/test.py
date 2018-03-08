from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def insert_read(session, model, file):
	with open(file, 'r') as file_fasta:
		i = 1
		sequence = ""
		for line in file_fasta:
			if ">" in line:
				sequence_id = line.strip()
				if i == 1:
					continue
				else:
					session.query(model).filter(model.sequence_id == sequence_id).update({model.read: sequence})
					sequence = ""
			else:
				sequence += line
			i += 1


def trim_reads(session, model):
	for element in session.query(model).all():
		read_sequence = element.read
		trimmed_read = element.read.replace(element.target, '')
		session.query(model).filter(model.sequence_id == element.sequence_id).update({model.read: trimmed_read})


def reverse_complement(session, model):
	for element in session.query(model).all():
		my_read = Seq(element.read, IUPAC.ambiguous_dna)
		reverse_read = my_read.reverse_complement()
		session.query(model).filter(model.sequence_id == element.sequence_id).update({model.reverse: reverse_read})



