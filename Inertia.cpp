#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

struct Cell
{
    int Index;
    vector<int> mon;
    vector<double> x;
    vector<double> y;
    vector<double> xCopy;
    vector<double> yCopy;
    double CenterX;
    double CenterY;
    vector<vector<double>> matrix;
    vector<double> eigenValue;
    vector<double> eigenVector;
    vector<double> eigenV;
    vector<double> Vp;
    int type;
    double eigenRatio;
    double area;
    double orderPcosa;

};

struct Monomer
{
    double x;
    double y;
    int type;

};
struct Area
{
    double Area1;
    double Area2;
};
struct ratioAnalysis
{
    vector<double> ratio1;
    vector<int> CellIndex1;
    vector<double> ratio2;
    vector<int> CellIndex2;
    double AveRatio1;
    double AveRatio2;

};

    vector<double> AR1;
    vector<double> AR2;

vector<Cell> Cells;
vector<Monomer> mon;
vector<Area> AveArea;

Cell Calinertia(Cell cell);
Cell CaleigenValue(Cell cell);
Cell CalRatio(Cell cell);
Cell CaleigenVector(Cell Cell);
Cell CalVposition(Cell Cell);
Area CalAveArea(vector<Cell> Cells);
double reverseBoundaryCoordinate(double x, double AreaLength);


int main(int argc, char* argv[])
{

    if(argc != 7)
    {
        cout << "usage: " << argv[0] << " Index CellPosition MonPosition FrameNumber DensityNum DensityCoefficient" <<endl;
        return 0;
    }

      double AreaLength = 400;
      vector<double> area;
      area.push_back(AreaLength);
      int densityNum = atoi(argv[5]);
      double DensityCoefficient = atof(argv[6]);
      int TotalLoop = atoi(argv[4]);
      int CellNum;
      int MonNum;
        if(TotalLoop%densityNum != 0)
{
    TotalLoop -= TotalLoop%densityNum;
}

int FrameperDensity = TotalLoop/densityNum;

vector<double> orderParameter;

      ifstream GetCellPos;
    GetCellPos.open(argv[2], ios::in);

    ifstream GetMonPos;
    GetMonPos.open(argv[3], ios::in);






 for(int lo = 0; lo < TotalLoop; lo++)
 {


    if(lo%FrameperDensity == 0&&lo != 0)
    {
        AreaLength *= DensityCoefficient;
        area.push_back(AreaLength);
    }
     //if(lo%1 == 0)
     if(lo != -1)
     {
    //cout << "Loop: " << lo << endl;
    char buf[1024];

    ratioAnalysis RatioAnalysis;





    /*
    for(int i = 0; i < Cells.size(); i++)
    {
        for(int j = 0; j < Cells[i].mon.size(); j++)
        {
            cout << Cells[i].mon[j] << " ";

        }
        cout << endl;
    }
*/

    //Read the position of Cells
    vector<double> temparr;
    double temppos;
    GetCellPos.getline(buf,1024);
    if(!atoi(buf))
    {
        GetCellPos.close();
        return 0;
    }
    CellNum = atoi(buf);
    Cells.resize(CellNum);
    for(int i = 0; i < CellNum; i++)
    {
        if(GetCellPos.eof())
        {
            cout << "File format not recognized!" << endl;
            GetCellPos.close();
            return 0;
        }
        GetCellPos.getline(buf,1024);
        stringstream input(buf);
        while(input >> temppos)
        {
            temparr.push_back(temppos);
        }
        Cells[i].Index = temparr[0];
        Cells[i].CenterX = temparr[1];
        Cells[i].CenterY = temparr[2];
        Cells[i].type = temparr[3];
        Cells[i].area = temparr[4];
        temparr.clear();
    }
    /*
    for(int i = 0; i < CellNum; i++)
    {
        cout << Cells[i].Index << " " << Cells[i].CenterX << " " << Cells[i].CenterY << " " << Cells[i].type << endl;
    }

*/

    //Read the position of Monomers
    GetMonPos.getline(buf,1024);
    if(!atoi(buf))
    {
        GetMonPos.close();
        return 0;
    }
    MonNum = atoi(buf);
    mon.resize(MonNum);
    GetMonPos.getline(buf,1024);
    for(int i = 0; i < MonNum; i++)
    {
        GetMonPos.getline(buf,1024);
        stringstream input(buf);
        while(input >> temppos)
        {
            temparr.push_back(temppos);
        }
        mon[i].type = temparr[0];
        mon[i].x = temparr[1];
        mon[i].y = temparr[2];
        temparr.clear();
    }
    /*
    for(int i = 0; i < mon.size(); i++)
    {
        cout << mon[i].type << " " << mon[i].x << " " << mon[i].y << endl;
    }
    cout << mon.size() << endl;
*/

    //Read the Index of monomer for each cell
    ifstream GetIndex;
    string InfileName = argv[1];
    stringstream getIname;
    getIname << InfileName;
    getIname >> InfileName;
    GetIndex.open(InfileName);
    int contCell = 0;
    int temp;
    while(!GetIndex.eof())
    {
        GetIndex.getline(buf,1024);
        stringstream input(buf);
        //cout << buf << endl;
        while(input >> temp)
        {
            Cells[contCell].mon.push_back(temp);

        }
        contCell++;
        //cout << contCell << endl;
    }







    //connect monomer and cell
    for(int i = 0; i < Cells.size(); i++)
    {
        for(int j = 0; j < Cells[i].mon.size(); j++)
        {
            //cout << "test1" << endl;
            Cells[i].x.push_back(mon[Cells[i].mon[j]].x);
            Cells[i].y.push_back(mon[Cells[i].mon[j]].y);
            //cout << "test2" << endl;
            //cout << Cells[i].x[j] << " " << Cells[i].y[j] << endl;
        }
    }


    //fix boundary condition
    for(int i = 0; i < Cells.size(); i++)
    {
        for(int j = 0; j < Cells[i].x.size(); j++)
        {
            if(fabs(Cells[i].CenterX - Cells[i].x[j]) < 0.5*AreaLength)
            {
                Cells[i].xCopy.push_back(Cells[i].x[j]);
            }else if(fabs(Cells[i].CenterX - Cells[i].x[j]) >= 0.5*AreaLength)
            {
                Cells[i].xCopy.push_back(reverseBoundaryCoordinate(Cells[i].x[j],AreaLength));
            }

            if(fabs(Cells[i].CenterY - Cells[i].y[j]) < 0.5*AreaLength)
            {
                Cells[i].yCopy.push_back(Cells[i].y[j]);
            }else if(fabs(Cells[i].CenterY - Cells[i].y[j]) >= 0.5*AreaLength)
            {
                Cells[i].yCopy.push_back(reverseBoundaryCoordinate(Cells[i].y[j],AreaLength));
            }

        }
    }

    /*
    cout << "start" << endl;
    for(int i = 0; i < Cells.size(); i++)
    {
        for(int j = 0; j < Cells[i].xCopy.size();j++)
        {
            cout << Cells[i].xCopy[j] << " " << Cells[i].yCopy[j] << endl;
        }
    }
    cout << " end" << endl;

*/


////////////////////////////////////////////////////////////////////////////////////////////

// Calculate inertia and ratio;
for(int i = 0; i < CellNum; i++)
{
    Cells[i] = Calinertia(Cells[i]);
    Cells[i] = CaleigenValue(Cells[i]);
    Cells[i] = CalRatio(Cells[i]);
    Cells[i] = CaleigenVector(Cells[i]);
    Cells[i] = CalVposition(Cells[i]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////
//define orientation, this part is for plot orientation graph
//if(lo%50 == 49 && lo != 499)
if(lo == -1)
{
ofstream oriout1,oriout2;
string orifname1,orifname2;
stringstream getoriName;
getoriName << "OriRig" << lo;
getoriName >> orifname1;
getoriName.clear();
getoriName << "OriSof" << lo;
getoriName >> orifname2;
getoriName.clear();
oriout1.open(orifname1,ios::app);
oriout2.open(orifname2,ios::app);



for(int i = 0; i < Cells.size();i++)
{
    //cout << Cells[i].Vp[0] << " " << Cells[i].Vp[1] << endl;
    //cout << Cells[i].Vp[2] << " " << Cells[i].Vp[3] << endl;
    if(Cells[i].type == 1)
    {
        oriout1 << Cells[i].Vp[0] << " " << Cells[i].Vp[1] << endl;
        oriout1 << Cells[i].Vp[2] << " " << Cells[i].Vp[3] << endl;
    }else if(Cells[i].type == 2)
    {
        oriout2 << Cells[i].Vp[0] << " " << Cells[i].Vp[1] << endl;
        oriout2 << Cells[i].Vp[2] << " " << Cells[i].Vp[3] << endl;
    }
}
}


////////////////////////////////////////////////////////
//output orientation for calculate relationship
ofstream outorientation;
string orienName;
stringstream oriFname;
oriFname << "orientation";
oriFname >> orienName;
outorientation.open(orienName, ios::app);
outorientation << Cells.size() << endl;
for(int i = 0; i < Cells.size();i++)
{
    outorientation << Cells[i].Index << " " << Cells[i].eigenV[1] << " " << Cells[i].eigenV[0] << endl; // Switch x and y at this time
}
outorientation.close();


//////////////////////////////////////////////////////////////////////////////////
//The part for calculate order parameter for each frame
//First, Find the average orientation for each frame
double AveOriX = 0, AveOriY = 0;
for(int i = 0; i < CellNum; i++)
{
    if(Cells[i].eigenV[0] < 0 && Cells[i].eigenV[1] < 0)
    {
        Cells[i].eigenV[0] = -Cells[i].eigenV[0];
        Cells[i].eigenV[1] = -Cells[i].eigenV[1];

    }else if(Cells[i].eigenV[0] < 0 && Cells[i].eigenV[1] >=  0)
    {
        Cells[i].eigenV[0] = -Cells[i].eigenV[0];
        Cells[i].eigenV[1] = -Cells[i].eigenV[1];
    }
    AveOriX += Cells[i].eigenV[0];
    AveOriY += Cells[i].eigenV[1];
}
double FrameVectorLength = sqrt(AveOriX*AveOriX + AveOriY*AveOriY);
orderParameter.push_back(0);
for(int i = 0; i < CellNum; i++)
{
    Cells[i].orderPcosa = (Cells[i].eigenV[0]*AveOriX + Cells[i].eigenV[1]*AveOriY)/(1*FrameVectorLength);
    Cells[i].orderPcosa = 2*Cells[i].orderPcosa*Cells[i].orderPcosa - 1;
    //cout << Cells[i].orderPcosa << endl;
    orderParameter[lo] += Cells[i].orderPcosa;
}
orderParameter[lo] /= CellNum;



/////////////////////////////////////////////////////////////////
//calculate average area
AveArea.push_back(CalAveArea(Cells));
/////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////
//Analysis result

if(lo != -1)
{
for(int i = 0; i < CellNum; i++)
{
    if(Cells[i].type == 1)
    {
        RatioAnalysis.ratio1.push_back(Cells[i].eigenRatio);
        RatioAnalysis.CellIndex1.push_back(i);
    }else if(Cells[i].type == 2)
    {
        RatioAnalysis.ratio2.push_back(Cells[i].eigenRatio);
        RatioAnalysis.CellIndex2.push_back(i);
    }
}

double avetemp = 0;
for(int i = 0; i < RatioAnalysis.ratio1.size(); i++)
{
    avetemp += RatioAnalysis.ratio1[i];
}
RatioAnalysis.AveRatio1 = avetemp;
RatioAnalysis.AveRatio1 /= RatioAnalysis.ratio1.size();

AR1.push_back(RatioAnalysis.AveRatio1);
avetemp = 0;
for(int i = 0; i < RatioAnalysis.ratio2.size(); i++)
{
    avetemp += RatioAnalysis.ratio2[i];
}
RatioAnalysis.AveRatio2 = avetemp;
RatioAnalysis.AveRatio2 /= RatioAnalysis.ratio2.size();

AR2.push_back(RatioAnalysis.AveRatio2);






/*

for(int i = 0; i < RatioAnalysis.ratio1.size();i++)
{
    cout << RatioAnalysis.ratio1[i] << endl;
}
    cout << "Averatio1: " << RatioAnalysis.AveRatio1 << endl;

for(int i = 0; i < RatioAnalysis.ratio2.size();i++)
{
    cout << RatioAnalysis.ratio2[i] << endl;
}
    cout << "Averatio2: " << RatioAnalysis.AveRatio2 << endl;


    for(int i = 0; i < Cells.size();i++)
    {
        cout << Cells[i].eigenRatio << endl;
        if(abs(Cells[i].eigenRatio - 1.01482) < 0.1)
        {
            cout << i << endl;
            cout << Cells[i].CenterX << " " << Cells[i].CenterY << endl;
            for(int j = 0; j < Cells[i].mon.size();j++)
            {
                cout << Cells[i].xCopy[j] << " " << Cells[i].yCopy[j] << endl;
            }
        }
    }


    // mark label
    cout << "Label1" << endl;
    if(lo !=0)
    {
    for(int i = 0; i < Cells.size();i++)
    {
        if(Cells[i].type == 1)
        {
        cout << Cells[i].CenterX << " " << Cells[i].CenterY << " " << Cells[i].eigenRatio << endl;
        //cout << "matrix" << endl;
        //cout << Cells[i].matrix[0][0] << " " << Cells[i].matrix[0][1] << endl;
        //cout << Cells[i].matrix[1][0] << " " << Cells[i].matrix[1][1] << endl;
        }
    }
    cout << "Label2" << endl;
    for(int i = 0; i < Cells.size();i++)
    {
        if(Cells[i].type == 2)
        {
        cout << Cells[i].CenterX << " " << Cells[i].CenterY << " " << Cells[i].eigenRatio << endl;
        //cout << "matrix" << endl;
        //cout << Cells[i].matrix[0][0] << " " << Cells[i].matrix[0][1] << endl;
        //cout << Cells[i].matrix[1][0] << " " << Cells[i].matrix[1][1] << endl;
        }
    }

    }
    */
     }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //do test
/*
    double Alpha = 2*3.1415926/80;
    Cell testring;
    testring.CenterX = 0;
    testring.CenterY = 0;
    for(int i = 0; i < 80; i++)
    {
        testring.mon.push_back(i);
        double rx = 5*cos(Alpha*i);
        double ry = 5*sin(Alpha*i);
        testring.xCopy.push_back(rx);
        testring.yCopy.push_back(ry);
        cout << testring.xCopy[i] << " " << testring.yCopy[i] << endl;
    }

    cout << endl;
    cout << "test" << endl;
    testring = Calinertia(testring);
    cout << "test" << endl;
    cout << testring.matrix[0][0] << " " << testring.matrix[0][1] << endl;
    cout << testring.matrix[1][0] << " " << testring.matrix[1][1] << endl;
    testring = CaleigenValue(testring);
    testring = CalRatio(testring);
    cout << testring.eigenRatio << endl;
*/
    /*
    if(testring.eigenValue[0] >= testring.eigenValue[1])
    {
        cout << testring.eigenValue[0]/testring.eigenValue[1] <<endl;
    }
    */


    Cells.clear();
    mon.clear();
 }
 }

 //Finish the loop


cout << "Inertial with time step"  << endl;
 for(int i = 0; i < AR1.size();i++)
 {
     cout << (AR1[i] + AR2[i])/2.0 << " " << AR1[i] << " " << AR2[i] << endl;
 }

/////////////////////////////////////////////////////////////
 //cal average Need change parameter
//Ave Inertia for many frames
 vector<double> AR1Ave;
 vector<double> AR2Ave;
 for(int i = 0; i < densityNum; i++)
 {
     AR1Ave.push_back(0);
     AR2Ave.push_back(0);
 }

 int aveCon = 0;
 for(int i = 0; i < AR1.size(); i++)
 {
     if(i % FrameperDensity == 0&&i != 0)
     {
         aveCon++;
     }
     AR1Ave[aveCon] += AR1[i];
     AR2Ave[aveCon] += AR2[i];
 }

 cout << endl;
 cout << "Average Inertia with density: " << endl;
 for(int i = 0; i < densityNum; i++)
 {
     AR1Ave[i] /= FrameperDensity;
     AR2Ave[i] /= FrameperDensity;
     //cout << 484/(area[i]*area[i]) << " " << AR1Ave[i] << " " << AR2Ave[i] << endl;
     cout << CellNum/(area[i]*area[i]) << " " << (AR1Ave[i]+AR2Ave[i])/2 << endl;
 }
//////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
// Area
cout << "Average Area with time:" << endl;
for(int i = 0; i < AveArea.size();i++)
{
    cout << i << " " << AveArea[i].Area1 << " " << AveArea[i].Area2 << endl;
}
int aveFrame = 0;
vector<Area> FrameArea;
FrameArea.resize(densityNum);
for(int i = 0; i < densityNum; i++)
{
    FrameArea[i].Area1 = 0;
    FrameArea[i].Area2 = 0;
}
int aveFcon = 0;
for(int i = 0; i < AveArea.size();i++)
{
    if(i % FrameperDensity == 0&&i != 0)
    {
        aveFcon++;
    }
    FrameArea[aveFcon].Area1 += AveArea[i].Area1;
    FrameArea[aveFcon].Area2 += AveArea[i].Area2;
}
cout << "Average Area for different density: " << endl;

    for(int i = 0; i < densityNum;i++)
    {
        FrameArea[i].Area1/= FrameperDensity;
        FrameArea[i].Area2/= FrameperDensity;
        //cout << 484/(area[i]*area[i]) << " " << FrameArea[i].Area1 << " " << FrameArea[i].Area2 << endl;
        cout << CellNum/(area[i]*area[i]) << " " << (FrameArea[i].Area1 + FrameArea[i].Area2)/2 << " " << endl;

    }





////////////////////////////////////////////////
//print order parameter
cout << " orderParameter VS time "<< endl;
for(int i = 0; i < orderParameter.size(); i++)
{
    cout << i << " " << orderParameter[i] << endl;
}

int orderPcount = 0;
vector<double> orderPvsDensity;
for(int i = 0; i < densityNum; i++)
{
    orderPvsDensity.push_back(0);
}

cout << " orderParameter Vs Density" << endl;
for(int i = 0; i < orderParameter.size(); i++)
{
    if(i %50 == 0 && i != 0)
    {
        orderPcount++;
    }
    orderPvsDensity[orderPcount] += orderParameter[i];
}
for(int i = 0; i < orderPvsDensity.size(); i++)
{
    orderPvsDensity[i] /= 50;
}

for (int i = 0; i < orderPvsDensity.size(); i++)
{
    cout << 484/(area[i]*area[i]) << " " << orderPvsDensity[i] << endl;
}

////////////////////////////////////////////////


    return 0;
}


Cell Calinertia(Cell cell)
{
    vector<vector<double>> temp;
    temp.resize(2);
    temp[0].resize(2);
    temp[1].resize(2);
    for(int i = 0; i < temp.size();i++)
    {
        for(int j = 0; j < temp[i].size(); j++)
        {
            temp[i][j] = 0;
        }
    }
    for(int i = 0; i < cell.mon.size(); i++)
        {

            temp[0][0] += pow((cell.yCopy[i] - cell.CenterY),2);
            temp[0][1] += (cell.xCopy[i] - cell.CenterX)*(cell.yCopy[i] - cell.CenterY);
            temp[1][1] += pow((cell.xCopy[i] - cell.CenterX),2);
        }
        temp[1][0] = temp[0][1];

        for(int i = 0; i <2; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                if(fabs(temp[i][j]) < 0.000001)
                {
                    temp[i][j] = 0;
                }
            }
        }

        cell.matrix = temp;



        return cell;

}

Cell CaleigenValue(Cell cell)
{
    double a, b, c;
    a = 1;
    b = -(cell.matrix[0][0] + cell.matrix[1][1]);
    c =  cell.matrix[0][0] * cell.matrix[1][1] - cell.matrix[0][1]*cell.matrix[1][0];
    vector<double> temp;
    double x1 = (-b + sqrt(pow(b,2) - 4*a*c))/(2*a);
    double x2 = (-b - sqrt(pow(b,2) - 4*a*c))/(2*a);

    temp.push_back(x1);
    temp.push_back(x2);
    cell.eigenValue = temp;
    return cell;
}



Cell CalRatio(Cell cell)
{
    cell.eigenRatio = (cell.eigenValue[0] > cell.eigenValue[1]) ? cell.eigenValue[0]/cell.eigenValue[1] : cell.eigenValue[1]/cell.eigenValue[0];
    return cell;
}


double reverseBoundaryCoordinate(double x, double AreaLength)
{
	if(x > 0.5*AreaLength)
	{
		x = x - AreaLength;
	}
	else if(x < 0.5*AreaLength)
	{
		x = x + AreaLength;
	}
	return x;
}


Cell CaleigenVector(Cell Cell)
{
    Cell.eigenVector.resize(2);
   double a,b,c,d,x1,x2;

        if(Cell.eigenValue[0] >= Cell.eigenValue[1])
        {
        a = Cell.matrix[0][0] - Cell.eigenValue[0];
        b = Cell.matrix[0][1];
        c = Cell.matrix[1][0];
        d = Cell.matrix[1][1] - Cell.eigenValue[0];
        }else if(Cell.eigenValue[0] < Cell.eigenValue[1])
        {
        a = Cell.matrix[0][0] - Cell.eigenValue[1];
        b = Cell.matrix[0][1];
        c = Cell.matrix[1][0];
        d = Cell.matrix[1][1] - Cell.eigenValue[1];
        }
        if(fabs(a) < 0.000001)
        {
            a = 0;
        }
        if(fabs(b) < 0.000001)
        {
            b = 0;
        }
        if(fabs(c) < 0.000001)
        {
            c = 0;
        }
        if(fabs(d) < 0.000001)
        {
            d = 0;
        }
        //cout << a << " " << b << endl;
        //cout << c << " " << d << endl;
        if(a!=0)
        {
            x2 = a/sqrt(a*a+b*b);
            x1 = -b/a*x2;
        }else if(c != 0)
        {
            x2 = c/sqrt(c*c + d*d);
            x1 = -d/c*x2;
        }else if(b != 0)
        {
            x1 = b/sqrt(a*a + b*b);
            x2 = -a/b*x1;
        }else if(d != 0)
        {
            x1 = d/sqrt(c*c + d*d);
            x2 = -c/d*x1;
        }

        if((a == 0&& d == 0&&c != 0&&d != 0)||(c == 0 && b == 0&&a!=0&&d!=0))
        {
            x1 = 0;
            x2 = 0;
        }
        if(a == 0 && b == 0 && c == 0 && d == 0)
        {
            x1 = 1;
            x2 = 0;
        }
        //cout << "x12" << endl;
        //cout << x1 << x2 << endl;
        //cout << "test2" << endl;
        Cell.eigenVector.push_back(x1);
        Cell.eigenVector.push_back(x2);

        //cout << "test3" << endl;
    double temp1,temp2;
    //temp1 = Cell.eigenVector[0][0] + Cell.eigenValue[1]/Cell.eigenValue[0]*Cell.eigenVector[1][0];
    //temp2 = Cell.eigenVector[0][1] + Cell.eigenValue[1]/Cell.eigenValue[0]*Cell.eigenVector[1][1];
    temp1 = x1;
    temp2 = x2;
    Cell.eigenV.push_back(temp1);
    Cell.eigenV.push_back(temp2);// x and y should switch in the real situation
    //cout << " pos" << endl;
    //cout << Cell.eigenV[0] << " " << Cell.eigenV[1] << endl;
    return Cell;
}


Cell CalVposition(Cell Cell)
{
    //cout << " pos" << endl;
    //cout << Cell.eigenV[0] << " " << Cell.eigenV[1] << endl;
    vector<double> temp;
    double a;
    a = Cell.CenterX - 0.1*0.5*Cell.eigenRatio*Cell.eigenV[1];
    temp.push_back(a);
    a = Cell.CenterY - 0.1*0.5*Cell.eigenRatio*Cell.eigenV[0];
    temp.push_back(a);
    a = Cell.CenterX + 0.1*0.5*Cell.eigenRatio*Cell.eigenV[1];
    temp.push_back(a);
    a = Cell.CenterY + 0.1*0.5*Cell.eigenRatio*Cell.eigenV[0];
    temp.push_back(a);

    Cell.Vp = temp;
    return Cell;
}

Area CalAveArea(vector<Cell> Cells)
{
    Area temp;
    temp.Area1 = 0;
    temp.Area2 = 0;
     int nth1 = 0;
     int nth2 = 0;
    for(int i = 0; i < Cells.size(); i++)
    {
        if(Cells[i].type == 1)
        {
            temp.Area1 += Cells[i].area;
            nth1++;
        }else if(Cells[i].type == 2)
        {
            temp.Area2 += Cells[i].area;
            nth2++;
        }
    }

    temp.Area1 /= nth1;
    temp.Area2 /= nth2;
    return temp;
}



