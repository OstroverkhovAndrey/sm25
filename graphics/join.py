
# Обединяет файлы полученные в результате работы MPI программы
# Пример вызова: $puthon3 join.py\
# ./../mpi/txt/M_80__N_120__MPI_2\ общий патерн файлов в которых лежат результаты работы каждого MPI процесса
#                                \ данные будут собираться их файлов M_80__N_120__MPI_2_0.txt  M_80__N_120__MPI_2_1.txt
# 2\ количество MPI процессов, оно же количество input файлов
# ./M_80__N_120__MPI_2.txt\ файл куда сохранять результат работы 

import argparse

def join(input, count, output):

    with open(output, "w") as outfile:
        for i in range(count):
            input_filename = input + f"_{i}.txt"
            with open(input_filename, "r") as infile:
                if i == 0:
                    for input_line in infile:
                        outfile.write(input_line)
                else:
                    for input_line in list(infile)[1:]:
                        outfile.write(input_line)
    print(f"Содержимое файлов успешно объединено в {output}")


def main():
    parser = argparse.ArgumentParser(description="Обединяет файлы полученные в результате работы MPI программы")
    parser.add_argument('input_filename_pattern', type=str, help='Паттерн назавний файлов от каждого MPI процесса')
    parser.add_argument('input_filename_count', type=int, help='Сколько всего файлов надо объединить')
    parser.add_argument('output_filename', type=str, help='Путь до файла куда сохранить реузльтат')
    args = parser.parse_args()

    join(args.input_filename_pattern, args.input_filename_count, args.output_filename)


if __name__ == '__main__':
    main()