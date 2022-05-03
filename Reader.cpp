#include <iostream>
#include <fstream> // To use ifstream
#include <vector>
using namespace std;

int main()
{
    std::vector<int> numbers;
    ifstream inputFile("Coarse.dat");        // Input file stream object

    // Check if exists and then open the file.
    if (inputFile.good()) {
        // Push items into a vector
        int current_number = 0;
        while (inputFile >> current_number){
            numbers.push_back(current_number);
        }

        // Close the file.
        inputFile.close();

        // Display the numbers read:
        cout << "The numbers are: ";
        for (int count = 0; count < numbers.size(); count++){
            cout << numbers[count] << " ";
        }

        cout << endl;
    }else {
        cout << "Error!";
        _Exit(0);
    }

    return 0;
}