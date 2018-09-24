//-----------------------------------------------------------------------------
// chan_map.c
// File to create channel map for BONuS12 RTPC
// Nate Dzbenski
// Last update: 12Jul2018
//-----------------------------------------------------------------------------

#include<math.h>
#include<vector>
#include<map>
#include<fstream>

using namespace std;


void chan_map()
{
    ofstream myfile;
    myfile.open ("chan_list.txt");
    
    // Establish constants
    const float PI = 3.1415926535;
    float PAD_W =2.79;
    float PAD_L = 4.0;
    float PAD_S = 80.0;
    float RTPC_L = 384.0;
    float phi_per_pad = PAD_W/PAD_S;
    
    int z_shift = 0;
    int Num_of_Col = RTPC_L/PAD_L;
    int Num_of_Row = (ceil) (2.0*PI*PAD_S)/PAD_W;
    int col = 0;
    int row = 0;
    int chan = 0;
    
    double z_pos = 0;
    double phi_rad =0;
    
    map< int, vector<double> > chanMap;
    vector<double> pos;
    
    chanMap.clear();
    
    const char separator    = '     ';
    const int width     = 12;
    
    myfile << left << setw(width) << setfill(separator) << "channel" << "z-pos[mm]" << left << setw(width) << setfill(separator) << "   phi-pos[rad]" << endl;
    
    // find row
    for(int i=0; i < Num_of_Row; i++){
        pos.clear();
        double phi_rad = (i*phi_per_pad)+(phi_per_pad/2.0);
        
        // find column
        for(int j=0; j < Num_of_Col; j++){
            row = (int) (phi_rad/phi_per_pad);
            float z_shift = row%4;
            
            z_pos = ((-RTPC_L/2.0 + PAD_L/2.0)+j*PAD_L) - z_shift;
            col = (int) (z_pos+z_shift+(RTPC_L/2.0))/PAD_L;
            chan = row*Num_of_Col+col;
            
            pos.push_back(z_pos);
            pos.push_back(phi_rad);
            
            // creates a output listing of z-position and phi-position for each channel
            myfile << left << setw(width) << setfill(separator) << chan << left << setw(width) << setfill(separator) << z_pos << left << setw(width) << setfill(separator) << phi_rad << endl;
            
            // creates a channel map with rows represent movement around phi [rad]
            // and columns represent movement in z [mm]
            //myfile << left << setw(width) << setfill(separator) << chan;
            
            chanMap.insert(make_pair (chan,pos));
            //cin.get();
        }
        //myfile << endl;
    }
    
    myfile.close();
    return;
}
