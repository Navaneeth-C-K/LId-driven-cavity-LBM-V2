//lid driven cavity flow
//L2 error crieterion has been implimented here!

//LID DRIVEN CAVITY FLOW SIMULATION -FINAL
//i am implimenting output to csv format!



//All the class stuff is going here !
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>

using namespace std;

class vector2D
{
    // double x, y;

    public:
        double x, y;
        vector2D()
        {
            x = 0;
            y = 0;
        }

        vector2D(double x_, double y_)
        {
            x = x_;
            y = y_;
        }

    double operator*(const vector2D &rhs)
    {
        return (this->x*rhs.x + this->y*rhs.y);
    }

    vector2D operator*(double val)
    {
        return vector2D(this->x*val, this->y*val);
    }

    vector2D operator+(const vector2D &rhs)
    {
        return vector2D(this->x + rhs.x, this->y + rhs.y);
    }

    vector2D operator=(const vector2D &rhs)
    {
        return vector2D(x = rhs.x, y = rhs.y);
    }

    vector2D operator/(double val)
    {
        if(val==0){
        cout << "Not defined" << endl;
        return vector2D(this->x, this->y);
        }
        else{
        return vector2D(this->x/val, this->y/val);
        }
    }
void print_vector()
    {
        cout << "X = " << x << ", Y = " << y <<endl;
    }

    friend ostream& operator<<(ostream &out, const vector2D &v);

    
};

ostream& operator<<(ostream &out, const vector2D &v)
    {
        out << "X = " << v.x << ", Y = " << v.y <<endl;

        return out;
    }

class direction
{
    //array<vector2D, 9> xi;

    public:
    array<vector2D, 9> xi;
        direction()
        {
            xi[0] = vector2D(0, 0);
            xi[1] = vector2D(1, 0);
            xi[2] = vector2D(0, 1);
            xi[3] = vector2D(-1, 0);
            xi[4] = vector2D(0, -1);
            xi[5] = vector2D(1, 1);
            xi[6] = vector2D(-1, 1);
            xi[7] = vector2D(-1, -1);
            xi[8] = vector2D(1, -1);
        }

        vector2D operator[](int rhs) const
        {
            return xi[rhs];
        }

};

class weight
{
    array<double, 9> w;

    public:
        weight()
        {
            w[0] = 4.0/9.0;
            w[1] = 1.0/9.0;
            w[2] = 1.0/9.0;
            w[3] = 1.0/9.0;
            w[4] = 1.0/9.0;
            w[5] = 1.0/36.0;
            w[6] = 1.0/36.0;
            w[7] = 1.0/36.0;
            w[8] = 1.0/36.0;


        }

        weight operator*(const double &rhs) 
        {
            weight product;
            for(int i = 0; i<9; i++)
            {
                product.w[i] =  w[i]*rhs;
            }
            return product;
        }

        double operator[](int rhs) const
        {
            return w[rhs];
        }
};

class node_property

{
    
    weight W;
    //direction Xi; will make this private later!
    
 

    public:
    direction Xi;
        array<double, 9> f;
        double rho;
        vector2D vel;

        

        void equilibrium(double rho, vector2D u, double Cs)//here i want to multiply w[i] with a variable rho 
        {

            for(int i= 0; i<9;i++)
            {
                f[i] =  (W[i] * rho)*(1 + ((u*Xi[i])/pow(Cs,2)) + (0.5* pow(u*Xi[i], 2)/pow(Cs,4)) - (0.5*(u*u)/pow(Cs, 2)) );
            }
        }

        double equilibrium(double rho, vector2D u, double Cs, int n) //this will help me create direction specific f_equilibrium (as it is impo for boundary node calculations)
        {
            //node_property equi; 

                return  (W[n] * rho)*(1 + ((u*Xi[n])/pow(Cs,2)) + (0.5* pow(u*Xi[n], 2)/pow(Cs,4)) - (0.5*(u*u)/pow(Cs, 2)) );
            
        }



};

//************************** */

//main code starts from here!

int main()
{

    //all the set up stuff happens here! 

    double Re = 400;
    int nodes = 529; //for Re = 400 i used 529 and Re = 100 i used 256
    double Ma = 0.1;
    double Cs = 1/sqrt(3.0);
    double U_star = Ma * Cs;
    double neu_star = (U_star * nodes) / Re;

    double Tau = (3*neu_star) + 0.5;

    double density = 1000;
    double rho = 1;
    double rho_constant = rho;
    vector2D vel_zero(0, 0);

    double total_time = 50.0;
    double L_x = 1.0;
    double U_lid = 1.0;

    int total_nodes = ((nodes + 2) * (nodes + 2));

    double node_length = L_x / nodes;
    double time_step = (U_star * node_length) / U_lid;
    


    cout << time_step << endl;
    cout << U_star << endl;

    vector<node_property> mesh(total_nodes);

    vector<node_property> temporary(total_nodes);

    vector<node_property> equilibrium_array(total_nodes);

    //setting equilibrium loop! this i think i have perfected!


    for(int i = 0; i < total_nodes; i++)
    {
        mesh[i].rho = rho;
        mesh[i].vel = vel_zero;
        mesh[i].equilibrium(mesh[i].rho, mesh[i].vel, Cs);

        //for(int j = 0; j < 9; j++)
        //{
          //  equilibrium_array[i].f[j] = mesh[i].f[j];
        //} this loop is completly useless as the equilibrium that we require the the f equation is dependent on the velocity at that particular time step
        
    } 


    int t = 0;
    double total_time_lu =  total_time/time_step;


    //***************************************************************************************** */

    //the time steps start from here! 

    vector<double>  u_array(total_nodes);
    vector<double>  v_array(total_nodes);

    double saved_squared_sum = 0;

    double error = 1;
    double error_u = 1;
    double error_v = 1;

    int steps = 0;

    while(error > 1e-6)
    {

        double squared_sum = 0.0;
        double squared_sum_v = 0.0;
        double sum_square_diff = 0.0;
        double sum_square_diff_v = 0.0;
        
        
        //**************************************************** */
        
        //the colission step! i think i have perfected this too (but this will still be something that i need to check out again but only the lowest priority!)

        for(int i = nodes+3; i <= total_nodes-(nodes+4); i++) //colission step
        {
            for(int j = 0; j < 9 ; j++)
            {
                if(i % (nodes + 2) == 0 || (i + 1) % (nodes + 2) == 0) {continue;}

                else
                {
                double N_d =  mesh[i].equilibrium(mesh[i].rho , mesh[i].vel, Cs, j);
                temporary[i].f[j] = mesh[i].f[j] - ((mesh[i].f[j]-N_d)/Tau); // here i am doing the collision and saving it as temporary array
                }

            }    
        }

    //********************************************* */
    //after we setup all the walls and corners we do the streaming ! this is also very very iffy i think this guy is the biggest source of my error!
    int n = 1;
    for(int i = nodes+3; i <= total_nodes-(nodes+4); i++) //streaming step
    {
        if(i % (nodes + 2) == 0 || (i + 1) % (nodes + 2) == 0) {continue;}

       else
       {
        for(int j = 0; j < 9 ; j++)
        {
            
            int a;
            if(j == 0){a = 0;}
            if(j == 1){a = 1;}
            if(j == 2){a = nodes+2;}
            if(j == 3){a = -1;}
            if(j == 4){a = -nodes-2;}
            if(j == 5){a = nodes +2 + 1;}
            if(j == 6){a = nodes + 2 - 1;}
            if(j == 7){a = -nodes - 2 - 1;}
            if(j == 8){a = -nodes - 2 + 1;}

            int target = i + a;

            mesh[target].f[j] = temporary[i].f[j];     
            //}
        }
        }
            
    }
    
        //************************************ */
        //i need a wall setup step and i should not have put streaming before the wall setup!.. so will put a wall set up loop!
        //i also made a horrible mistake that is i put streaming inside the wall setup which is WRONG!

        //Boundary set up!
        
        for(int i = nodes+3; i <= total_nodes-(nodes+4); i++)
        {
            if(i % (nodes + 2) == 0 || (i + 1) % (nodes + 2) == 0) {continue;}

            else if(i ==  total_nodes-(nodes+nodes+3)) //top left corner
            {
                double u_x_top = U_star;
                double u_y_top = 0;
                double u_x_left = 0;
                double u_y_left = 0; 
            
                mesh[i].rho = (mesh[i+1].rho + mesh[i-nodes].rho)/2;

                vector2D corner_vel = vector2D((u_x_left + u_x_top)/2, (u_y_left + u_y_top)/2); 
            
                mesh[i].f[4] = mesh[i].f[2] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 2) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 4); 
            
                mesh[i].f[1] = mesh[i].f[3] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 3) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 1);  
            
                mesh[i].f[5] = (mesh[i].rho*(1 + (u_y_left + u_y_top)/2) - mesh[i].f[0] - mesh[i].f[1] - (2* mesh[i].f[2]) - mesh[i].f[3] - (2*mesh[i].f[6]))*0.5;

                mesh[i].f[7] = ((-mesh[i].rho*(((u_x_left + u_x_top)/2) + ((u_y_left + u_y_top)/2))) + mesh[i].f[1] + ( mesh[i].f[2]) - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]) )*0.5;

                mesh[i].f[8] = (-mesh[i].rho*((u_y_left + u_y_top)/2)) + mesh[i].f[2] - mesh[i].f[4] + mesh[i].f[5] + mesh[i].f[6] - mesh[i].f[7];


            }

            else if(i == total_nodes-(nodes+4)) //top right corner
            {

                double u_x_top = U_star;
                double u_y_top = 0;
                double u_x_right = 0;
                double u_y_right = 0; 

                mesh[i].rho = (mesh[i-1].rho + mesh[i-nodes].rho)/2;

                vector2D corner_vel = vector2D((u_x_right + u_x_top)/2, (u_y_right + u_y_top)/2);

     
                mesh[i].f[3] = mesh[i].f[1] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 1) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 3);
            
                mesh[i].f[4] = mesh[i].f[2] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 2) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 4); 
 
                mesh[i].f[8] = ((mesh[i].rho*(1+(u_x_right+u_x_top)/2)) - mesh[i].f[0] - (2*mesh[i].f[1]) - mesh[i].f[2] -mesh[i].f[4] - (2*mesh[i].f[5]) )*0.5;

                mesh[i].f[7] = ((-mesh[i].rho*((u_x_right + u_x_top)/2) + ((u_y_right + u_y_top)/2)) + mesh[i].f[1] + ( mesh[i].f[2]) - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]) )*0.5; 

                mesh[i].f[6] = (mesh[i].rho*(u_y_right + u_y_top)) - mesh[i].f[2] + mesh[i].f[4] - mesh[i].f[5] + mesh[i].f[7] + mesh[i].f[8];
            

            }

            else if(i == nodes+3) //bottom left corner
            {


                double u_x_bottom = 0;
                double u_y_bottom = 0;
                double u_x_left = 0;
                double u_y_left = 0;

                mesh[i].rho = (mesh[nodes].rho + mesh[1].rho)/2;

                vector2D corner_vel = vector2D((u_x_left + u_x_bottom)/2, (u_y_left + u_y_bottom)/2);

                mesh[i].f[1] = mesh[i].f[3] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 3) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 1);

                mesh[i].f[2] = mesh[i].f[4] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 4) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 2);

                mesh[i].f[8] = (mesh[i].rho*(1 - (u_y_left + u_y_bottom)/2) - mesh[i].f[0] - mesh[i].f[1] - mesh[i].f[3] - (2*mesh[i].f[4]) - (2*mesh[i].f[7]))*0.5;

                mesh[i].f[5] = ((mesh[i].rho*((u_x_left + u_x_bottom)/2) + ((u_y_left + u_y_bottom)/2)) - mesh[i].f[1] - (mesh[i].f[2]) + mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[7]) )*0.5; 
            
                mesh[i].f[6] = (mesh[i].rho*( (u_y_left + u_y_bottom)/2)) - mesh[i].f[2] + mesh[i].f[4] - mesh[i].f[5] + mesh[i].f[8] + mesh[i].f[7];
            
            

            }

            else if(i == nodes+2+nodes) //bottom right corner
            {

            double u_x_bottom = 0;
            double u_y_bottom = 0;
            double u_x_right = 0;
            double u_y_right = 0;

            mesh[i].rho = (mesh[nodes-2].rho + mesh[nodes + nodes-1].rho)/2;

            vector2D corner_vel = vector2D((u_x_right + u_x_bottom)/2, (u_y_right + u_y_bottom)/2);

            mesh[i].f[3] = mesh[i].f[1] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 1) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 3);

            mesh[i].f[2] = mesh[i].f[4] - mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 4) + mesh[i].equilibrium(mesh[i].rho, corner_vel, Cs, 2);

            mesh[i].f[7] = (mesh[i].rho*(1 - (u_y_right + u_y_bottom)/2) - mesh[i].f[0] - mesh[i].f[1] - mesh[i].f[3] - (2*mesh[i].f[4]) - (2*mesh[i].f[8]))*0.5;

            mesh[i].f[6] = ((mesh[i].rho*(1 - (u_x_right + u_x_bottom)/2)) - mesh[i].f[0] - mesh[i].f[2] - (2*mesh[i].f[3]) - mesh[i].f[4] - (2*mesh[i].f[7]))*0.5; //had forgotten to multiply this with 0.5 !

            mesh[i].f[5] = (mesh[i].rho*((u_y_right + u_y_bottom)/2)) - mesh[i].f[2] + mesh[i].f[4] - mesh[i].f[6] + mesh[i].f[7] + mesh[i].f[8];
            

            }

            else if((i+2)%(nodes+2) == 0)//right wall
            {
            int a = 1;
            int b = 0;

            double u_x = 0;
            double u_y = 0;

            vector2D wall_vel(u_x, u_y); 

            double f_equal = mesh[i].f[0] + mesh[i].f[2] + mesh[i].f[4];
            double f_plus = mesh[i].f[5] + mesh[i].f[1] + mesh[i].f[8];
            
            mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));

            mesh[i].f[3] = mesh[i].f[1] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 1) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs,  3); //here when i am trying to find f_eq that is wrong i think! as the i am using Vel_zero but as we know the top wall does have some velocity 

            mesh[i].f[6] = ((mesh[i].rho*(u_y - u_x)) + mesh[i].f[1] - mesh[i].f[2] - mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[8]))*0.5;
            mesh[i].f[7] = ((mesh[i].rho*(-u_x-u_y)) + mesh[i].f[1] + mesh[i].f[2] - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]))*0.5;

            
            
            }

            else if((i-1)%(nodes+2) == 0)//left wall
            {
            int a = -1;
            int b = 0;

            double u_x = 0;
            double u_y = 0;

            vector2D wall_vel(u_x, u_y); 

            double f_equal = mesh[i].f[0] + mesh[i].f[2] + mesh[i].f[4];
            double f_plus = mesh[i].f[6] + mesh[i].f[3] + mesh[i].f[7];
            
            mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));

            mesh[i].f[1] = mesh[i].f[3] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 3) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 1);

            mesh[i].f[5] = ((mesh[i].rho*(u_x + u_y)) - mesh[i].f[1] - mesh[i].f[2] + mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[7]))*0.5;
            mesh[i].f[8] = ((mesh[i].rho*(u_x - u_y)) - mesh[i].f[1] + mesh[i].f[2] + mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[6]))*0.5;


            
            }

            else if(i<(total_nodes-(nodes + 4)) && i >(total_nodes-(nodes + nodes + 3)))//top wall
            {
            int a = 0;
            int b = 1;

            double u_x = U_star; //some thing to solve!!!!!!
            double u_y = 0;

            vector2D wall_vel(u_x, u_y); 
             
            double f_equal = mesh[i].f[0] + mesh[i].f[1] + mesh[i].f[3];
            double f_plus = mesh[i].f[6] + mesh[i].f[2] + mesh[i].f[5];
            
            mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));

            mesh[i].f[4] = mesh[i].f[2] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 2) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs,  4);

            mesh[i].f[7] = ((mesh[i].rho*(-u_x - u_y)) + mesh[i].f[1] + mesh[i].f[2] - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]))*0.5;
            mesh[i].f[8] = ((mesh[i].rho*(u_x - u_y)) - mesh[i].f[1] + mesh[i].f[2] + mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[6]))*0.5;
            


            }

            else if(i>(nodes+3) && i<(nodes+2+nodes))//bottom wall
            {


            int a = 0;
            int b = -1;

            double u_x = 0;
            double u_y = 0;

            vector2D wall_vel(u_x, u_y); 

            double f_equal = mesh[i].f[3] + mesh[i].f[0] + mesh[i].f[1];
            double f_plus = mesh[i].f[7] + mesh[i].f[4] + mesh[i].f[8];
            
            mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));

            mesh[i].f[2] = mesh[i].f[4] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 4) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs,  2);

            mesh[i].f[5] = ((mesh[i].rho*(u_x + u_y)) - mesh[i].f[1] - mesh[i].f[2] + mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[7]))*0.5;
            mesh[i].f[6] = ((mesh[i].rho*(-u_x + u_y)) + mesh[i].f[1] - mesh[i].f[2] - mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[8]))*0.5;

            
            }

  

        }
    

        //***************************** */ */
        //and finally the loop to calculate all the rho's and velocity vectors this i think i have perfected but still on the caution list!
        for(int i = nodes+3; i <= total_nodes-(nodes+4); i++)
        {
            if(i % (nodes + 2) == 0 || (i + 1) % (nodes + 2) == 0) {continue;} 
            else{
            vector2D sum_f_xi_temp(0, 0); 
            
            
            mesh[i].rho = 0;
            for(int j = 0; j < 9 ; j++)//this loop will just be there for velocity calculation!
            {
                mesh[i].rho = mesh[i].rho + mesh[i].f[j]; // i initially set the special loop just for rho calculation here but i guess i could just calculate the rho when streaming!
                    
                
                sum_f_xi_temp = sum_f_xi_temp + (mesh[i].Xi[j] * mesh[i].f[j]);
            }

            
            mesh[i].vel = sum_f_xi_temp/mesh[i].rho;

            if(steps > 1000 && steps % 100 == 0)
            {    
                sum_square_diff = sum_square_diff + pow((mesh[i].vel.x - u_array[i]), 2);
                sum_square_diff_v = sum_square_diff_v + pow((mesh[i].vel.y - v_array[i]), 2);

                squared_sum = squared_sum + pow((mesh[i].vel.x), 2);
                squared_sum_v = squared_sum_v + pow((mesh[i].vel.y), 2);
            }

            u_array[i] = mesh[i].vel.x;
            v_array[i] = mesh[i].vel.y;
        }
            
        }

        cout << mesh[nodes+3].rho << endl;
        cout << mesh[nodes+3+nodes].rho << endl;
        cout << mesh[total_nodes-nodes-3].rho <<endl;
        cout << mesh[total_nodes -nodes-nodes-2].rho <<endl;
        cout <<mesh[35].rho <<endl;

        //i am just commenting these are these will take up precious procession power for nothing...it was useful for finding errors...but other than that these are quite useless!
        // double sum=0;
        // for(int i = nodes+3; i < total_nodes-(nodes+4); i++)
        // {
        //     sum = sum + mesh[i].rho;
        // }
        // cout <<"*********"<<endl;
        // cout << sum << endl;
        cout<<"runnig!" <<endl;

        cout << t << "/" << total_time_lu << endl;
        if(steps > 1000 && steps % 100 == 0)
        {
            // error_u = sqrt(sum_square_diff / squared_sum); these are rmemenants of the past! i was thinking of doing error_u + error_v / 2 but that is not proper L2 stuff!!
            // error_v = sqrt(sum_square_diff_v / squared_sum_v);

            error = sqrt((sum_square_diff + sum_square_diff_v)/(squared_sum + squared_sum_v));

            cout << "Error = " << error <<endl;
        }

        t++;
        steps++;

        saved_squared_sum  = squared_sum;
        saved_squared_sum = squared_sum_v;
    }
for(int i = nodes+3; i < total_nodes-(nodes+4); i++)
{
cout << mesh[i].rho << endl;
}
vector<vector2D> vel(total_nodes);
for(int i= 0; i<total_nodes; i++)
{
cout << mesh[i].rho << endl;
vel[i] = mesh[i].vel*U_star;
}

// for(int i =  total_nodes-(nodes+nodes+3); i <= total_nodes-(nodes+4); i++)
// {
//     cout << mesh[i].vel << endl;
// }



int pad = nodes + 2; //i am just adding the two graveyard nodes on both sides...now when dividing we get the center perfectly 
int center_x = pad / 2; 

ofstream outfile1("centerline_y.csv");

outfile1 << "y_norm,u_x,u_y\n";

for(int i = 0; i < total_nodes; i++)
{
    int x = i % pad;
    int y = i / pad;

    
    if(x == 0 || x == pad - 1 || y == 0 || y == pad - 1) { continue; }

    vector2D temp_vel;
    if(x == center_x)
    {

        double y_norm = (double)(y - 1) / (nodes - 1);

        temp_vel = mesh[i].vel/U_star ;
        outfile1 << y_norm << "," << temp_vel.x << "," << temp_vel.y << "\n";
    }

    
}
    outfile1.close();
    cout << "Data successfully written to centerline_y.csv!" << endl;

    cout << "****************************************************************** " << endl;

        //int pad = nodes + 2; //i am just adding the two graveyard nodes on both sides...now when dividing we get the center perfectly 
        int center_y = pad / 2; 

    
        int start_i = (center_y * pad) + 1;
        int end_i = start_i + nodes - 1;

        int j = 0;
        vector2D temp_vel;

        ofstream outfile2("centerline_x.csv");

        outfile2 << "x_norm,u_x,u_y\n";
    
        for(int i = start_i; i <= end_i; i++)   
        {
            double x_norm = (double)j / (nodes - 1);
            temp_vel = mesh[i].vel / U_star;
        
            outfile2 << x_norm << "," << temp_vel.x << "," << temp_vel.y << "\n";
            j++;
        }

        outfile2.close();
        cout << "Data successfully written to centerline_x.csv!" << endl;

//----------------------------------------------------------------------------------------------------------------------
//**************************TESTING */

    ofstream outField("velocity_field.csv");

    outField << "x,y,u,v\n";

   

    // Loop only over the actual fluid domain (1 to nodes)

    for(int y = 1; y <= nodes; y++) {

        for(int x = 1; x <= nodes; x++) {

            // Calculate the 1D index from the 2D coordinates

            int i = (y * pad) + x;

            outField << x << "," << y << "," << mesh[i].vel.x/U_star << "," << mesh[i].vel.y/U_star << "\n";
            //outField << x << "," << y << "," << mesh[i].vel.x << "," << mesh[i].vel.y << "density =" << mesh[i].rho << "\n";


        }

    }

    outField.close();

    cout << "2D velocity field saved to velocity_field.csv!" << endl;

}