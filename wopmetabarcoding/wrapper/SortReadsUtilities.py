import os
from wopmars.utils.Logger import Logger
from wopmetabarcoding.utils.utilities import get_or_create
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import sqlite3


def create_primer_tag_fasta_for_vsearch(sample_information_obj, forward_strand, primer_tag_fasta):
    """
    Function creating a fasta file which will be used as a database by vsearch

    :param session: Current session of the database
    :param sample_information_obj: Model of the SampleInformation table
    :param primer_tag_fasta:
    :param forward_strand: Boolean
    :return:
    """
    with open(primer_tag_fasta, 'w') as fout:
        for row in sample_information_obj:
            # if row.tag_forward != "" and row.primer_forward != "" and row.tag_reverse != "" and row.primer_reverse != "":
            if forward_strand:
                fout.write(">" + row.tag_forward + row.primer_forward + "\n")
                fout.write(row.tag_forward + row.primer_forward + "\n")
            else:
                fout.write(">" + row.tag_reverse + row.primer_reverse + "\n")
                fout.write(row.tag_reverse + row.primer_reverse + "\n")


# def read_counter(session, file, model):
#     """
#     Function counting occurences of a read in the fasta file and store it in a table
#     :param session: Current of the database
#     :param file: fasta containing the reads
#     :param model: Model of the ReadCount table
#     :return: void
#     """
#     with open(file, "r") as fasta_file:
#         next(fasta_file)
#         sequence = ""
#         liste_tmp = []
#         for line in fasta_file:
#             if ">" in line:
#                 if sequence in liste_tmp:
#                     session.query(model).filter(model.sequence == sequence).update({model.count: model.count+1})
#                     sequence = ""
#                 else:
#                     obj_readcount = {"sequence": sequence, "count": 1}
#                     insert_table(session, model, obj_readcount)
#                     liste_tmp.append(sequence)
#                     sequence = ""
#             else:
#                 sequence += line.strip()
#         if session.query(model.sequence.contains(sequence)) is True and sequence != "":
#             session.query(model).filter(model.sequence == sequence).update({model.count: model.count + 1})
#             sequence = ""
#         else:
#             obj_readcount = {"sequence": sequence, "count": 1}
#             insert_table(session, model, obj_readcount)
#             sequence = ""


def check_criteria_in_vsearch_output(vsearch_output_tsv, checked_vsearch_output_tsv, overhang_value):
    """
    Function to select lines of vsearch_output_tsv according to these SRS criteria:
    1. first position of target in the alignment (tilo) is 1
    2. last position of target in the alignment (tihi) equals length of the target (tl)
    3. lower case part of the target (tag) has an exact match (non-case sensitive) to the aligned part of query (qrow)
    4. this match is within -overhang bases in the 5’ end of qrow

    Example of vsearch_output_csv:
    M00842:118:000000000-ABGKE:1:1101:19978:2032|tgatcgatgatcag|AGATATTGGAACWTTATATTTTATTTTTGG	atatcgacagtgWACTAATCAATTWCCAAATCCTCC	36	1	36	1	36	ATATCGACAGTGTACTAATCAATTTCCAAATCCTCC
    M00842:118:000000000-ABGKE:1:1101:19826:1953|acgatccacagtg|AGATATTGGAACWTTATATTTTATTTTTGG	tctcgatgatcagWACTAATCAATTWCCAAATCCTCC	37	1	37	1	37	TCTCGATGATCAGAACTAATCAATTTCCAAATCCTCC

    :param session: Current session of the database
    :param model: Model of the table
    :param vsearch_output_tsv: csv containing the results of the alignment
    :return: void
    """
    with open(checked_vsearch_output_tsv, 'w') as fout:
        with open(vsearch_output_tsv, 'r') as fin:
            next(fin)
            i = 0
            for line in fin:
                if line.split("\t")[5] == "1":
                    tag_primer_sequence_query = line.split('\t')[1]
                    vsearch_tl = line.split('\t')[2]
                    vsearch_tilo = line.split('\t')[5]
                    vsearch_tihi = line.split('\t')[6]
                    tag_primer2read_alignment = line.split('\t')[7].strip()
                    tag_sequence = ""
                    for nucleotide in tag_primer_sequence_query:
                        if nucleotide.islower():
                            tag_sequence += nucleotide.upper()
                    tag_position = tag_primer2read_alignment.find(tag_sequence)
                    if vsearch_tilo == "1" and vsearch_tihi == vsearch_tl \
                            and tag_sequence in tag_primer2read_alignment:
                        fout.write(line)
                    else:
                        i += 1
            Logger.instance().info(str(i) + " reads were discarded.")


def fasta_writer(session, model, file_name):
    """
    Function
    :param session:
    :param model:
    :param file_name:
    :return:
    """
    name = file_name.replace(".fasta", "_count.fasta")
    with open(name, "w") as fout:
        j = 1
        line = []
        for element in session.query(model).all():
            fout.write("> Read number " + str(j) + " count: " + str(element.count))
            fout.write("\n")
            fout.write(element.sequence + '\n')
            j += 1


def trim_reads(checked_vsearch_output_tsv, read_fasta_file_name, trimmed_out_tsv, is_forward_strand, tempdir):
    """

    :param checked_vsearch_output_tsv:
    :param read_fasta_file_name:
    :param session:
    :param file_model:
    :param is_forward_strand:
    :return:
    """
    #
    # Import read_fasta_file_name into sqlite for indexing and faster access
    # read_fasta_db_sqlite = read_fasta_file_name.replace('.fasta', '.sqlite')
    read_fasta_db_sqlite = os.path.join(tempdir, os.path.basename(read_fasta_file_name).replace('.fasta', '.sqlite'))
    conn = sqlite3.connect(read_fasta_db_sqlite)
    try:
        conn.execute("DROP TABLE IF EXISTS read_fasta")
        conn.execute("CREATE TABLE  read_fasta (id VARCHAR PRIMARY KEY , seq VARCHAR)")
        for record in SeqIO.parse(read_fasta_file_name, 'fasta'):
            conn.execute("INSERT INTO read_fasta (id, seq) VALUES (?, ?)", (str(record.description.split()[0]), str(record.seq)))
        conn.commit()
    except UnicodeDecodeError:
        pass
    with open(trimmed_out_tsv, 'w') as fout:
        with open(checked_vsearch_output_tsv, 'r') as fin:
            for line in fin:
                line_info = line.strip().split('\t')
                read_cursor = conn.execute('SELECT seq FROM read_fasta WHERE id=?', (line_info[0],))
                read_list = list(read_cursor.fetchone())
                read_cursor.close()
                read = read_list[0]
                vsearch_qihi = line_info[4]
                primer_tag_seq = read[0:int(vsearch_qihi)]
                trimmed_read = read.replace(primer_tag_seq, "")
                trimmed_read_seq_obj = Seq(trimmed_read, IUPAC.ambiguous_dna)
                reverse_read = trimmed_read_seq_obj.reverse_complement()
                line = line.strip() + '\t' + str(reverse_read) + '\n'
                fout.write(line)
    conn.close()


def convert_trimmed_tsv_to_fasta(trimmed_tsv, trimmed_fasta):
    # new_fasta = forward_trimmed_fasta.replace('.csv', '.fasta')
    # csv_file = open(forward_trimmed_fasta, 'r')
    with open(trimmed_tsv, 'r') as fin:
        with open(trimmed_fasta, 'w') as fout:
            for line in fin:
                # line = line.strip()
                line = line.strip().split("\t")
                tag = "".join([character for character in line[1] if character.islower()])
                marker = "".join([character for character in line[1] if character.isupper()])
                fout.write(">" + line[0] + "|" + tag + "|" + marker + "\n")
                fout.write(line[8] + "\n")
    # csv_file.close()


def annotate_reads(session, file_information_model, trimmed_tsv, file_id, annotated_reads_tsv):
    """
    Function used to merge all file information between the vsearch results reads and the file information csv
    :param session: Current session of the database
    :param file_information_model: SampleInformation table, contains all the necessary information like marker_id names
    :param model2: File table, contains all the information about files used
    :param trimmed_tsv: csv file to parse
    :param merged_fasta_file_name: Name of the original merged fasta file
    :return: void
    """
    # # Creating a merged_fasta_file_name for the output csv file
    # annotated_reads_tsv = trimmed_tsv.replace('_forward_trimmed_output_reverse_trimmed.csv', '_combination.tsv')
    # # Inserting the merged_fasta_file_name in the database for used it later in an iteration
    # session.query(model2).filter(model2.file_name == merged_fasta_file_name).update({model2.final_csv: annotated_reads_tsv})
    # # Open output file
    # output = open(annotated_reads_tsv, 'w')
    with open(annotated_reads_tsv, 'w') as fout:
        # Parsing input csv file
        with open(trimmed_tsv, 'r') as fin:
            i = 0
            incoherent_list = []
            for line in fin:
                line = line.strip().split('\t')
                #
                id_tag_primer = line[0]
                read_id = id_tag_primer.split('|')[0]
                tag_forward = id_tag_primer.split('|')[1]
                tag_reverse = "".join([character for character in line[1] if character.islower()])
                read_sequence = line[8]
                # Get the information in the SampleInformation table with tags and helped by the merged_fasta_file_name
                file_information_instance = session.query(file_information_model)\
                    .filter(file_information_model.tag_forward == tag_forward)\
                    .filter(file_information_model.tag_reverse == tag_reverse)\
                    .filter(file_information_model.file_id == file_id).first()
                if file_information_instance is not None:

                    try:
                        # Write a new tsv with these informations
                        fout.write(
                            read_id + "\t"
                            + str(file_information_instance.marker_id) + "\t"
                            + file_information_instance.run_name +"\t"
                            + file_information_instance.tag_forward + "\t"
                            + file_information_instance.tag_reverse + "\t"
                            + file_information_instance.sample_name + "\t"
                            + file_information_instance.replicate_name + "\t"
                            + read_sequence + "\n"
                        )
                    except AttributeError:
                        print(file_information_instance.tag_forward, file_information_instance.tag_reverse, file_information_instance.marker_id)
                else:
                    i += 1
                    incoherent_list.append(read_id)
            Logger.instance().info(str(i) + " incoherent reads discarded.")


def count_reads(gathered_marker_file, count_reads_marker, marker_name, sample_count_tsv):
    """

    :param gathered_marker_file:
    :param count_reads_marker:
    :return:
    """
    # Parse the database for Marker files
    # Creating a database name to store the results
    # Store the database file filename in the marker_model for later use
    # Open connection with the database file
    conn = sqlite3.connect(count_reads_marker)
    # Drop the table if the program has been launch before
    conn.execute("DROP TABLE IF EXISTS count_read")
    # Create a table in the databse file
    conn.execute("CREATE TABLE  count_read (id VARCHAR , marker_id INT, count INT, seq VARCHAR PRIMARY KEY)")
    biosamples_count_variant = {}
    biosample_count = {}
    # Parse the csv
    with open(gathered_marker_file, 'r') as marker_file:
        i = 1
        for line in marker_file:
            line = line.strip()
            line = line.split('\t')
            if line[7] in biosamples_count_variant:
                temp_count_per_samplereplicate = biosamples_count_variant.get(line[7])
                sample_replicate = line[5] + "_" + line[6]
                if sample_replicate in temp_count_per_samplereplicate:
                    temp_count_per_samplereplicate[sample_replicate] += 1
                else:
                    temp_count_per_samplereplicate[sample_replicate] = 1
                biosamples_count_variant[line[7]] = temp_count_per_samplereplicate
            else:
                sequence = line[7]
                sample_replicate = line[5] + "_" + line[6]
                count_per_samplereplicate = biosample_count.copy()
                count_per_samplereplicate[sample_replicate] = 1
                biosamples_count_variant[sequence] = count_per_samplereplicate
            check_read = conn.execute('SELECT EXISTS (SELECT id FROM count_read WHERE seq=?)', (line[7],))
            for row in check_read.fetchone():
                #  In case of the line already exist, the count number is updated
                if row != 0:
                    conn.execute('UPDATE count_read SET count = count + 1 WHERE seq=?', (line[7],))
                    # Else a row is created
                else:
                    variant_id = line[1] + "_variant_" + str(i)
                    marker_name = line[1]
                    read_count = 1
                    sequence = line[7]
                    conn.execute("INSERT INTO count_read (id, marker_id, count, seq) VALUES (?, ?, ?, ?)", (variant_id, marker_name, int(read_count), sequence))
                    i += 1
            check_read.close()
        # Line which delete singletons
        with open(sample_count_tsv, 'w') as test_output:
            test_output.write("{}\t{}\t{}\t{}\t{}\n".format('sequence', 'replicate', 'sample','sample_replicate', 'count'))
            for variant_seq in biosamples_count_variant:
                count_per_samplerepl = biosamples_count_variant.get(str(variant_seq))
                for sample_repl in count_per_samplerepl:
                    count = count_per_samplerepl.get(str(sample_repl))
                    if len(count_per_samplerepl) == 1 and count == 1:
                        continue
                    else:
                        replicate = sample_repl.split('_')[-1]
                        sample = sample_repl.split('_')[-2]
                        test_output.write(variant_seq + '\t' + replicate + '\t'+ sample + '\t' + sample_repl + '\t' + str(count))
                        test_output.write('\n')
        conn.execute("DELETE FROM count_read WHERE count=1")
        conn.commit()
        conn.close()


def gather_files(marker_name, gathered_marker_file, annotated_tsv_list, run_list):
    """
    Function used to gather all files of the same marker_id in one to count the reads occurrences
    :param marker_name: Marker name
    :param gathered_marker_file: Output file
    :param annotated_tsv_list: List of annotated files
    :return: void
    """
    # Search the csv files for each marker_id in the Marker table
    # for element in session.query(marker_model).all():
    # Creating a filename for a database file for each marker_id
    # filename = element.name + "_file"
    # file_place = "data/" + filename + ".csv"
    # Storing this name in the row of the corresponding marker_id
    # session.query(marker_model).filter(marker_model.name == element.name).update({marker_model.marker_file: file_place})
    i = 0
    for annotated_file in annotated_tsv_list:
        if i == 0:
            previous_run = run_list.get(annotated_file)
            gathered_marker_file = gathered_marker_file.replace(".tsv", "_" + previous_run + ".tsv")
            file = open(gathered_marker_file, 'w')
            file.close()
        if marker_name in annotated_file and previous_run in annotated_file:
            with open(gathered_marker_file, 'a') as fout:
                with open(annotated_file, 'r') as input_file:
                    fout.write(input_file.read())
        previous_run = run_list.get(annotated_file)
        i += 1


def insert_variant(session, count_reads_marker, variant_model):
    """
    Function used to insert variants in the variant table
    :param session: Current session of the database
    :param count_reads_marker: Temp database with a table for count database
    :param variant_model: Variant table
    :return: void
    """
    # Opening the database
    conn = sqlite3.connect(count_reads_marker)
    # Selecting some attributes in the database files
    variant_data = conn.execute("SELECT id, marker_id, seq, count FROM count_read")
    for row in variant_data.fetchall():
        obj_variant = {'variant_id': row[0], 'marker_id': row[1], 'sequence': row[2], 'readcount': row[3]}
        # Inserting theses information in the variant table
        instance = session.query(variant_model).filter(variant_model.sequence == row[2]).first()
        if not instance:
            get_or_create(session, variant_model, **obj_variant)
    variant_data.close()
    conn.close()











