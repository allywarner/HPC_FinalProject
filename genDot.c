
//Generates the dot file to view in graphviz

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lanczos.h"

int coordCompare (const void* a, const void* b);

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
    fscanf(readFile,"%*d %*d %*d\n");

    //reads and writes until the end of the read file, writes a color based on the processor
    int node1,node2,process;
    while(fscanf(readFile,"%d %d %d\n",&node1,&node2,&process) != EOF) {
        if(node1 == node2){
            continue;
        }
        fprintf(dotFile,"%d -- %d [style=filled,fillcolor = %s,fixedsize=true,color=%s];\n",node1,node2,colors[process],colors[process]);
    }

    //Writes the last line
    fprintf(dotFile,"}");

    //Closes the files
    fclose(readFile);
    fclose(dotFile);
}

void coord2Dot(coord* connections, int coordLength, int process){

    FILE* dotFile;

    //Opens new file to write
    dotFile = fopen("dotFile.gc","a");
    int i;

    for(i = 0;i < coordLength;i++){
        if(connections[i].row < connections[i].col){
            int temp = connections[i].row;
            connections[i].row = connections[i].col;
            connections[i].col = temp;
        }
    }

    //sort coordinates by row
    quicksort(connections,coordLength,sizeof(coord),coordCompare);

    for(i=0;i<coordLength;i+=2)
      printf("%d %d %d\n", connections[i].row,connections[i].col,process);

    //Writes the file for the length of the struct, writes a color based on the processor
    int node1,node2;
    for (i = 0; i < coordLength; i+=2) {
        node1 = connections[i].row;
        node2 = connections[i].col;
        fprintf(dotFile,"%d -- %d [style=filled,fillcolor = %s,fixedsize=true,color=%s];\n",node1,node2,colors[process],colors[process]);
    }

    //Closes the file
    fclose(dotFile);

}

int coordCompare (const void* a, const void* b) {
  struct point {
    int x;
    int y;
  };

  if ((*(struct point*)a).x > (*(struct point*)b).x)
    return 1;
  if ((*(struct point*)a).x == (*(struct point*)b).x) {
    if ((*(struct point*)a).y > (*(struct point*)b).y)
      return 1;
    else if ((*(struct point*)a).y == (*(struct point*)b).y)
      return 0;
    else
      return -1;
  }
  else /*(a < b)*/
    return -1;
}

// ////This main is for debugging
// int main(int argc, char* argv[]) {
//
//     if (argc < 2){
//         printf("Error.");
//         return 1;
//     }
//
//     char* file = argv[1];
//     file2Dot(file);
// }
