
   //
   // Sagan Simulator (Alpha Version) - a gravitational simulator code.
   // by Filipo Novo Mór - filipo.mor at gmail.com
   // 2011, November
   // this still a premature and test code yet!
   //

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TRUE  (-1)
#define FALSE (0)

// reserved for future use.
typedef struct
{
    long Code;
    char Description[50];
    char Source[20];
} error_type;
error_type GlobalError;

typedef struct
{
    long    ID;
    double  Mass;
    float   PosX;
    float   PosY;
    float   PosZ;
    float   vX;
    float   vY;
    float   vZ;
    double  fX;
    double  fY;
    double  fZ;
    char    Processed; // reserved for future use.
} particuleData;


typedef struct
{
    particuleData Content;
     struct particule* Next;
     struct particule* Prev;
     //struct particule* Parent;
} particule;


int FileExists (char *FileName)
{
	FILE *file_p;

    printf("... trying to check the file %s...\n", FileName);
	file_p = fopen(FileName, "r");
	if(file_p)
	{
		fclose(file_p);
		printf("... file '%s' found!!!!...\n", FileName);
		return(TRUE);
	}
	return(FALSE);
}

particule LineParsing(char *LineBuffer)
{
	particule Element;

	int  Cont      = 0;
	int  PosAux    = 0;
	int  Position  = 0;
	int  LineSize;
	char BufferAux[500];

	LineSize = strlen(LineBuffer);

	printf("==> Line being processed: %s\n", LineBuffer);

	// parse the line, finding the fields separated by blank spaces " " until reach the end of the line.
	while (Position <= LineSize)
	{

		if(LineBuffer[Position] != ' ') // consume the line, character by character.
		{
			BufferAux[Position-PosAux] = LineBuffer[Position];
		}
		if((LineBuffer[Position] == ' ' ) || (Position==LineSize)) // field separator or end of the line???
		{
			BufferAux[Position-PosAux] = NULL;
			switch (++Cont)
			{
				case 1: Element.Content.PosX = atof(BufferAux); break;
				case 2: Element.Content.PosY = atof(BufferAux); break;
				case 3: Element.Content.PosZ = atof(BufferAux); break;
				case 4: Element.Content.vX   = atof(BufferAux); break;
				case 5: Element.Content.vY   = atof(BufferAux); break;
				case 6: Element.Content.vZ   = atof(BufferAux); break;
				case 7: Element.Content.Mass = atof(BufferAux); Cont=0; break;
			}
			PosAux = Position+1;
		}
		Position++;
	}
	return Element;
}

void InsertNew(particule *ParentNode,
              long  ID,
              float PosX,
              float PosY,
              float PosZ,
              float vX,
              float vY,
              float vZ,
              double fX,
              double fY,
              double fZ,
              float Mass,
              char  Processed )
{
        // create a new node
        particule *NewNode;
        NewNode = (particule*)malloc(sizeof(particule));

        NewNode->Content.ID        = ID;
        NewNode->Content.PosX      = PosX;
        NewNode->Content.PosY      = PosY;
        NewNode->Content.PosZ      = PosZ;
        NewNode->Content.vX        = vX;
        NewNode->Content.vY        = vY;
        NewNode->Content.vZ        = vZ;
        NewNode->Content.fX        = fX;
        NewNode->Content.fY        = fY;
        NewNode->Content.fZ        = fZ;
        NewNode->Content.Mass      = Mass;
        NewNode->Content.Processed = Processed;
        NewNode->Prev = NULL;
        NewNode->Next = NULL;



        // list is empty.
        if(ParentNode->Next == NULL)
        {
            ParentNode->Next = NewNode;
            NewNode->Prev    = ParentNode;
            NewNode->Next    = NULL;
        }
        // list has at least one node
        else
        {
            // insert the node node immediately under the parent.
            // it's more efficient than run through over the whole list, right?
            particule *AuxNode;

            AuxNode = ParentNode->Next;
            NewNode->Next = AuxNode;
            AuxNode->Prev = NewNode;
            NewNode->Prev = ParentNode;
            ParentNode->Next = NewNode;
        }
}

void PrintList(particule *ParentNode)
{
    particule *Node;

    printf("\n>>>Listing nodes...\n");
    for(Node = ParentNode->Next; Node != NULL; Node = Node->Next)
    {
        printf("ID:   %d\n", Node->Content.ID);
        printf("****************\n");
        printf("PosX: %f\n",   Node->Content.PosX);
        printf("PosY: %f\n",   Node->Content.PosY);
        printf("PosZ: %f\n",   Node->Content.PosZ);
        printf("vX:   %f\n",   Node->Content.vX);
        printf("vY:   %f\n",   Node->Content.vY);
        printf("vZ:   %f\n",   Node->Content.vZ);
        printf("Mass: %f\n",   Node->Content.Mass);
        if(Node->Content.Processed==TRUE)
        {
              printf("Proc: TRUE\n\n");
        }else printf("Proc: FALSE\n\n");

    }
}

void Log(char* Buffer)
{

    //
    // this is a function used for test purpouses only
    //

    FILE *file_p;

    file_p = fopen("log.txt", "a");

    fprintf(file_p, Buffer);

    fclose(file_p);
}


void SaveTestFile(char* FileName, particule *ParentNode)
{

    //
    // this is a function used for test purpouses only
    //

    FILE *file_p;
    particule *Node;

    file_p = fopen(FileName, "w");

    fprintf(file_p, "\n>>>Listing nodes...\n");

    for(Node = ParentNode->Next; Node != NULL; Node = Node->Next)
    {
        fprintf(file_p, "ID:   %d\n", Node->Content.ID);
        fprintf(file_p, "****************\n");
        fprintf(file_p, "PosX: %f\n",   Node->Content.PosX);
        fprintf(file_p, "PosY: %f\n",   Node->Content.PosY);
        fprintf(file_p, "PosZ: %f\n",   Node->Content.PosZ);
        fprintf(file_p, "vX:   %f\n",   Node->Content.vX);
        fprintf(file_p, "vY:   %f\n",   Node->Content.vY);
        fprintf(file_p, "vZ:   %f\n",   Node->Content.vZ);
        fprintf(file_p, "fX:   %f\n",   Node->Content.fX);
        fprintf(file_p, "fY:   %f\n",   Node->Content.fY);
        fprintf(file_p, "fZ:   %f\n",   Node->Content.fZ);
        fprintf(file_p, "Mass: %f\n",   Node->Content.Mass);
        if(Node->Content.Processed==TRUE)
        {
              fprintf(file_p, "Proc: TRUE\n\n");
        }else fprintf(file_p, "Proc: FALSE\n\n");

        printf("\n --> saving particle %d...\n", Node->Content.ID);

    }
    fclose(file_p);
}

void GenerateResultFile(char* FileName, particule *ParentNode)
{

    //
    // used to generated the data file containing the updated particles positions.
    //

    FILE *file_p;
    particule *Node;

    file_p = fopen(FileName, "w");

    for(Node = ParentNode->Next; Node != NULL; Node = Node->Next)
    {
        fprintf(file_p, "%f %f %f %f %f %f %f\n", Node->Content.PosX, Node->Content.PosY, Node->Content.PosZ, 
		                                          Node->Content.vX, Node->Content.vY, Node->Content.vZ, 
												  Node->Content.Mass);
		//fprintf(file_p, "%f %f %f\n", Node->Content.PosX, Node->Content.PosY, Node->Content.PosZ);
        //fprintf(file_p, "%d %d %d\n", Node->Content.PosX, Node->Content.PosY, Node->Content.PosZ);
    }

    fclose(file_p);

}

int LoadFile(char *FileName, particule *ListHead)
{

    //
    // Load the file containing the initial conditions.
    //

	FILE *file_p;
	char  LineBuffer[100];
	int   Contador;

	printf("\Opening file %s... \n", FileName);

	if(!FileExists(FileName))
	{
		printf("\n File %s doesn't exists!!! \n", FileName);
		return(FALSE);
	}

	file_p = fopen(FileName, "r");
	if(!file_p)
	{
		printf("\nError opening data file !\n");
		return(FALSE);
	}

	Contador=0;
	while( fgets(LineBuffer, 100, file_p) != NULL )
	{
        particule P;
        P = LineParsing(LineBuffer);
        P.Content.ID = Contador++;
        P.Next = NULL;
        P.Prev = NULL;

        InsertNew(ListHead,
                  P.Content.ID,
                  P.Content.PosX,
                  P.Content.PosY,
                  P.Content.PosZ,
                  P.Content.vX,
                  P.Content.vY,
                  P.Content.vZ,
                  0.0, /* force on X axis */
                  0.0, /* force on Y axis */
                  0.0, /* force on Z axis */
                  P.Content.Mass,
                  P.Content.Processed);

	};

	close(file_p);

	printf("\n %s was closed successfully!\n", FileName);

	return(TRUE);
};

double cube(double value)
{
    return (value * value * value);
}

double sqr(double value)
{
    return ( value * value);
};

void CalculateForce(particule* PartBase, particule* PartNext)
{

    //
    // Calculate the force between two particules (PartBase and PartNext)
    //

    double Force;
    double XDist;
    double YDist;
    double ZDist;
    double XT;
    double YT;
    double ZT;
    double DistanceSQR;
    double Rinv;
    double Rinv3;

    // WARNING!!! Don't forget to update this line with the amount of particules to be
    // processed! eps = 4/n  !!!! Use decimal point! eps is the softening factor!
    //double eps = 4.0/395.0;
	double eps = 0;

    XDist = PartBase->Content.PosX - PartNext->Content.PosX;
    XT = sqr(XDist);

    YDist = PartBase->Content.PosY - PartNext->Content.PosY;
    YT = sqr(YDist);

    ZDist = PartBase->Content.PosZ - PartNext->Content.PosZ;
    ZT = sqr(ZDist);

    DistanceSQR = XT + YT + ZT;

    Rinv = 1.0 / sqrt(DistanceSQR + sqr(eps));

    Rinv3 = cube(Rinv);


    PartBase->Content.fX -= PartNext->Content.Mass * XDist * Rinv3;
    PartBase->Content.fY -= PartNext->Content.Mass * YDist * Rinv3;
    PartBase->Content.fZ -= PartNext->Content.Mass * ZDist * Rinv3;

}

double CalculaNCorpos(particule* ParentNode, double dt)
{

    //
    //  The particles are saved in a linked list, generated from the load of the initial conditions file.
    //  This function scans the whole linked list, calculating the gravitational attraction between each particle
    // with each other, one by one. This way, we will have N^2 interactions.
    //  After comparing all particles, finally update the positions and speeds.
    //  It's a time stepped algorithm, so each linked list scan represents a step in time, whose evolution value is
    // given by the 'dt' parameter. It's MANDATORY to keep this value fewer than '1'.
    //

    particule* NodeI;
    particule* NodeJ;

    for(NodeI = ParentNode->Next; NodeI != NULL; NodeI = NodeI->Next)
    {
        // updating the coordinates
        NodeI->Content.PosX += (NodeI->Content.vX * dt)/2;
        NodeI->Content.PosY += (NodeI->Content.vY * dt)/2;
        NodeI->Content.PosZ += (NodeI->Content.vZ * dt)/2;
    }


    for(NodeI = ParentNode->Next; NodeI != NULL; NodeI = NodeI->Next)
    {
        NodeI->Content.fX = 0.0;
        NodeI->Content.fY = 0.0;
        NodeI->Content.fZ = 0.0;

        for(NodeJ = ParentNode->Next; NodeJ != NULL; NodeJ = NodeJ->Next)
            if(NodeI->Content.ID != NodeJ->Content.ID) // we don't have to calculate the node against itself, right?
            {
                    CalculateForce(NodeI, NodeJ);
            }
    }


    for(NodeI = ParentNode->Next; NodeI != NULL; NodeI = NodeI->Next)
    {

        // updating the velocities
        NodeI->Content.vX   += NodeI->Content.fX * dt;
        NodeI->Content.vY   += NodeI->Content.fY * dt;
        NodeI->Content.vZ   += NodeI->Content.fZ * dt;

        // updating the coordinates
        NodeI->Content.PosX += (NodeI->Content.vX * dt)/2;
        NodeI->Content.PosY += (NodeI->Content.vY * dt)/2;
        NodeI->Content.PosZ += (NodeI->Content.vZ * dt)/2;

    }
    return dt;

}

int main()
{

    // create and set up the linked list head node.
    particule *ListHead;
    ListHead = malloc(sizeof(particule));
    ListHead->Prev = NULL;
    ListHead->Next = NULL;

    // WARNING! Don't forget to update here the name of the initial conditions file!
	char InitialConditionsFileName[] = "earth_jupiter_sun.ini";

	char ResultFileName[100];
	char NroArq[] = "  ";

    LoadFile(InitialConditionsFileName, ListHead);


    int    i = 0;          // we need these sequencial numbers to name the result files.
    double t = 0.0;        // controls the amount of time steps already processed for the simulation.
    double dt = 0.0001;    // indicates the contribution of each scan for the simulation evolution.
    double dtout = 0.01;   // indicates the value of the step where a result file has to be generated.
    double tout = t+dtout;
    int    timeSteps = 10; //WARNING!! Don't forget to update this value with the desired time steps amount!

    while (t<=timeSteps)
    {

        t += CalculaNCorpos(ListHead, dt);
        if (t>tout)
        {
            tout = t + dtout;

            sprintf(ResultFileName, "%03d.dat", i++);
            GenerateResultFile(ResultFileName, ListHead);
        }

    }

    return 0;
}
