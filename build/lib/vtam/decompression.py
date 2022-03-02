from utils.FileDecompression import FileDecompression

file = "../../merge/fastagz/mfzr_2_fw.fasta.gz"

decompress = FileDecompression(file)
fileDecompressed = decompress.gzip_decompression()
print(fileDecompressed)
decompress.delete_file()

