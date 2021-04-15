#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <chrono>
#include <vector>
#include <cmath>

//globals
const int size = 7000000;       // size % 2 = 0
const int Tag = 0;
const int p = 2;

void show_array(int* array, int size_ar);

int comp(const int* i, const int* j);

int* copy_array(int* array);



int main(int argc, char* argv[])
{

    int rc;
    int i, j, k;

    if (rc = MPI_Init(&argc, &argv))
    {
        std::cout << "Ошибка запуска, выполнение остановлено " << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    int rank;
    int numprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    //int array[size] = { 2,3,3,8,5,6,1,4 };
    //int* array = new int[size + (numprocs - size % numprocs)] { -222211111,-222211110,-22221119,-22221118,-22221117,-22221116,-22221115,-22221114,-22221113,-22221112,-22221111, -22221111,-22221110,-2222119,-2222118,-2222117,-2222116,-2222115,-2222114,-2222113,-2222112,-2222111, -2222111,-2222110,-222219,-222218,-222217,-222216,-222215,-222214,-222213,-222212,-222211, -222211,-222210,-22229,-22228,-22227,-22226,-22225,-22224,-22223,-22222,-22221, -11111,-11110,-1119,-1118,-1117,-1116,-1115,-1114,-1113,-1112,-1111, -1111,-1110,-119,-118,-117,-116,-115,-114,-113,-112,-111, -111,-110,-19,-18,-17,-16,-15,-14,-13,-12,-11, -11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1 };
    int* array = new int[size + (numprocs - size % numprocs)];
    for (i = 0; i < size + (numprocs - size % numprocs); i++)
        array[i] = std::rand();

    /*
    if (rank == 0)
    {
        std::cout << "Initial size: " << size << std::endl;
        std::cout << "% size: " << size % numprocs << std::endl;
        std::cout << "Initial array: " << std::endl;
        show_array(array, size + (numprocs - size % numprocs));
    }
    */

    MPI_Status status;

    auto start = std::chrono::high_resolution_clock::now();

    int num_proc_arr = (size + (numprocs - size % numprocs)) / numprocs;

    int* proc_arr = new int[num_proc_arr * 2]{};
    int* temp_proc_arr = new int[num_proc_arr]{};

    /*
    for (i = 0; rank + i * numprocs < (size + (numprocs - size % numprocs)); i++)
    {
        proc_arr[i] = array[rank + i * numprocs];
        if (rank == p)
            std::cout << "SSS size: " << proc_arr[i] << std::endl;
    }
    */

    for (i = 0; i < num_proc_arr && i + rank * num_proc_arr < (size + (numprocs - size % numprocs)); i++)
        proc_arr[i] = array[i + rank * num_proc_arr];




    /*
    if (size / numprocs == num_proc_arr)
    {
        int* proc_arr = new int[num_proc_arr * 2];
        int* temp_proc_arr = new int[num_proc_arr];

        for (i = 0; i < num_proc_arr; i++)
            proc_arr[i] = array[i + rank * num_proc_arr];
    }
    else
    {
        int* proc_arr = new int[(num_proc_arr + 1)* 2];
        int* temp_proc_arr = new int[num_proc_arr + 1];

        for (i = 0; i + rank * (num_proc_arr + 1) < size; i++)
            proc_arr[i] = array[i + rank * (num_proc_arr + 1)];

        if (rank == numprocs - 1)
            show_array(proc_arr, num_proc_arr + 1);
    }*/




    for (i = 1; i <= numprocs; i++)
    {
        if (i % 2 != 0)
        {
            for (j = 0; j + 1 < numprocs; j+=2)
            {
                if (rank == j + 1)
                    MPI_Send(proc_arr, num_proc_arr, MPI_INT, j, Tag, MPI_COMM_WORLD);

                if (rank == j)
                {
                    MPI_Recv(temp_proc_arr, num_proc_arr, MPI_INT, j + 1, Tag, MPI_COMM_WORLD, &status);

                    for (k = num_proc_arr; k < 2 * num_proc_arr; k++)
                        proc_arr[k] = temp_proc_arr[k - num_proc_arr];

                    qsort(proc_arr, 2 * num_proc_arr, sizeof(int), (int(*) (const void*, const void*)) comp);



                    for (k = num_proc_arr; k < 2 * num_proc_arr; k++)
                    {
                        temp_proc_arr[k - num_proc_arr] = proc_arr[k];
                        proc_arr[k] = 0;
                    }

                    MPI_Send(temp_proc_arr, num_proc_arr, MPI_INT, j + 1, Tag, MPI_COMM_WORLD);
                }

                if (rank == j + 1)
                {
                    MPI_Recv(temp_proc_arr, num_proc_arr, MPI_INT, j, Tag, MPI_COMM_WORLD, &status);

                    for (k = 0; k < num_proc_arr; k++)
                        proc_arr[k] = temp_proc_arr[k];
                }
            }
        }
        else
        {
            for (j = 1; j + 1 < numprocs; j += 2)
            {
                if (rank == j + 1)
                    MPI_Send(proc_arr, num_proc_arr, MPI_INT, j, Tag, MPI_COMM_WORLD);

                if (rank == j)
                {
                    MPI_Recv(temp_proc_arr, num_proc_arr, MPI_INT, j + 1, Tag, MPI_COMM_WORLD, &status);

                    for (k = num_proc_arr; k < 2 * num_proc_arr; k++)
                        proc_arr[k] = temp_proc_arr[k - num_proc_arr];

                    qsort(proc_arr, 2 * num_proc_arr, sizeof(int), (int(*) (const void*, const void*)) comp);



                    for (k = num_proc_arr; k < 2 * num_proc_arr; k++)
                    {
                        temp_proc_arr[k - num_proc_arr] = proc_arr[k];
                        proc_arr[k] = 0;
                    }

                    MPI_Send(temp_proc_arr, num_proc_arr, MPI_INT, j + 1, Tag, MPI_COMM_WORLD);
                }

                if (rank == j + 1)
                {
                    MPI_Recv(temp_proc_arr, num_proc_arr, MPI_INT, j, Tag, MPI_COMM_WORLD, &status);

                    for (k = 0; k < num_proc_arr; k++)
                        proc_arr[k] = temp_proc_arr[k];
                }
            }
        }
    }

    if (rank != 0)
        MPI_Send(proc_arr, num_proc_arr, MPI_INT, 0, Tag, MPI_COMM_WORLD);

    if (rank == 0)
    {
        int* new_array = new int[size + (numprocs - size % numprocs)];

        for (j = 0; j < num_proc_arr; j++)
        {
            new_array[j] = proc_arr[j];
        }

        for (i = 1; i < numprocs; i++)
        {
            MPI_Recv(temp_proc_arr, num_proc_arr, MPI_INT, i, Tag, MPI_COMM_WORLD, &status);
            for (j = 0; j < num_proc_arr; j++)
            {
                new_array[j + i * num_proc_arr] = temp_proc_arr[j];
            }
        }
        //show_array(new_array, size + (numprocs - size % numprocs));
    }

    auto end_first = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration_first = end_first - start;

    if (rank == 0)
        std::cout << "Duration Shell  sort = " << duration_first.count() << "sec" << std::endl;

    MPI_Finalize();

    /*
    //Bubble
    auto start = std::chrono::high_resolution_clock::now();

    int* bubble_array = bubble_sort(array, argc, argv);

    auto end_first = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration_first = end_first - start;

    std::cout << std::endl;
    std::cout << "Bubble sort: " << std::endl;
    show_array(bubble_array, size);
    std::cout << "Duration Bubble sort = " << duration_first.count() << "sec" << std::endl;
    */
    return 0;
}



int comp(const int* i, const int* j)
{
    return *i - *j;
}

void show_array(int* array, int size_ar)
{
    for (int i = 0; i < size_ar; i++)
    {
        std::cout.width(5);
        std::cout << array[i];
    }
    std::cout << std::endl;
}

int* copy_array(int* array)
{
    int* copy = new int[size];
    for (int i = 0; i < size; i++)
        copy[i] = array[i];
    return copy;
}
