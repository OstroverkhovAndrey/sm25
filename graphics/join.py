
output_filename = "all.txt"

with open(output_filename, "w") as outfile:
    for i in range(32):
        filename = f"./txt/32_{i}.txt"
        with open(filename, "r") as infile:
            outfile.write(infile.read())
            # Если нужно разделять содержимое файлов, можно добавить перевод строки:
            # outfile.write("\n")
print(f"Содержимое файлов успешно объединено в {output_filename}")