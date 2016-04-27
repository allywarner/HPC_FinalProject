//Generates the dot file to view in graphviz

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct _coord{
    int row;
    int col;
}coord;

//Initializes 16 colors. These can be changed to anything and I tried to pick colors that look different enough! More colors can be found at http://www.graphviz.org/doc/info/colors.html if you don't like them!
char* colors[] = {"black","red","blue","purple","yellow","magenta","brown","green","orange","skyblue1","forestgreen","slateblue1","wheat4","lightpink2","ivory4","cyan"};

//Generates the dot file with either one color or no color, we can easily change this if we are getting a file with several processors worth of information.
void file2Dot(char* data){

    //Initializes files
    FILE* dotFile;
    FILE* readFile;

    //Opens the file to read and write
    readFile = fopen(data,"r");
    dotFile = fopen("dotFile.gv","w");

    //Writes the first two lines
    fprintf(dotFile,"graph {\n");
    fprintf(dotFile,"node [shape = point]\n");

    //Gets rid of the first line in the file we are reading because it does not have connection or processor data
    fscanf(readFile,"%*[^\n]\n",NULL);

    //reads and writes until the end of the read file, writes a color based on the processor
    int node1,node2,process;
    while(fscanf(readFile,"%d %d %d\n",&node1,&node2,&process) != EOF) {
        fprintf(dotFile,"%d -- %d [style=filled,fillcolor = %s,fixedsize=true,color=%s];\n",node1,node2,colors[process],colors[process]);
    }

    //Writes the last line
    fprintf(dotFile,"}");

    //Closes the files
    fclose(readFile);
    fclose(dotFile);
}

void coord2Dot(coord* connections, size_t n, int process){

    //Initializes file
    FILE* dotFile;

    //Opens new file to write
    dotFile = fopen("dotFile.gc","w");

    //Writes the first two lines
    fprintf(dotFile,"graph {\n");
    fprintf(dotFile,"node [shape = point]\n");

    //Writes the file for the length of the struct, writes a color based on the processor
    int node1,node2;
    for (int i = 0; i < n; i++) {
        node1 = connections[i].row;
        node2 = connections[i].col;
        fprintf(dotFile,"%d -- %d [style=filled,fillcolor = %s,fixedsize=true,color=%s];\n",node1,node2,colors[process],colors[process]);
    }

    //Writes the last line
    fprintf(dotFile,"}");

    //Closes the file
    fclose(dotFile);

}

////This main is for debugging
int main(int argc, char* argv[]) {

    if (argc < 2){
        printf("Error.");
        return 1;
    }

    char* file = argv[1];
    file2Dot(file);
}
