//Generates the dot file to view in graphviz

#include "lanczos.h"

char* colors[] = {"green","red","blue","purple","yellow","pink","brown"};

void genDot(char* data,int colorNum){
    
    FILE* dotFile,readFile;
    
    readFile = fopen(data,"r");
    dotFile = fopen("dotFile.gv","w");
    
    fprintf(dotFile,"graph {");
    
    int node1,node2;
    while(fscanf(readFile,"%d %d\n",&node1,&node2) != EOF) {
        if (colorNum == 0) {
        fprintf(dotFile,"%d -- %d;\n",node1,node2);
        } else {
            fprintf(dotFile,"%d -- %d [color = %s];",node1,node2,colors[colorNum-1])
        }
    }
    fprintf(dotFile,"}");
    
    fclose(readFile);
    fclose(dotFile);
}

//This main is for debugging
int main(int argc, char* argv[]) {
    
    if (argc < 2){
        cerr << "Error.";
        return 1;
    }
    
    char* file = argv[1];
    int color = atoi(argv[2]);
    genDot(file,color);
}