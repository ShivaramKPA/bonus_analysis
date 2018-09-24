import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

// usage:
// rungroovy Eloss.groovy <filename.hipo>

HipoDataSource reader = new HipoDataSource();

// define plots


// open file
reader.open(args[0]);

// check for events
while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
        if(event.hasBank("RTPC::pos" && "RTPC::adc")){
        
        // read in banks
        DataBank bank_adc = event.getBank("RTPC::adc");
        DataBank bank_pos = event.getBank("RTPC::pos");
        
        // read in variables from banks
        
        
        
        
        // calculations
        
        
        }
        
        
        
        
        
        
}
        


// 

