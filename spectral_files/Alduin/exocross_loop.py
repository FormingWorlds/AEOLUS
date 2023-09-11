import os
import subprocess

def xcross_compute_xsec(temperature_range, pressure_range, wavenumber_range, spectral_resolution, clean):
    
    path_data     = "/proj/bolinc/users/x_ryabo/spectral_files/ExoMol/raw_data/"
    path_exocross = "/proj/bolinc/users/x_ryabo/spectral_files/ExoMol/exocross/" 
    path_inputs   = "/proj/bolinc/users/x_ryabo/spectral_files/ExoMol/input_files/"
    path_outputs  = "/proj/bolinc/users/x_ryabo/spectral_files/ExoMol/xsec/"
    if not os.path.exists(path_inputs):
        os.makedirs(path_inputs)
    if not os.path.exists(path_outputs):
        os.makedirs(path_outputs)
        
    # Loop over temperature and pressure
    for temperature in temperature_range:
        for pressure in pressure_range:
            # Create a unique input file name based on temperature and pressure
            input_file_name  = path_inputs+f"H2O_{temperature}K_P{pressure}bar.inp"
            output_file_name = path_outputs+f"1H2-16O_{int(temperature):04d}K_Voigt_P{pressure:.3e}_cstep_"+'{:03d}'.format(int(float(spectral_resolution) * 100))+"hires.xsec"

            # Create the input file
            with open(input_file_name, "w") as input_file:
                input_file.write(f"Temperature {temperature:.1f} (K)\n")
                input_file.write(f"Range 0.0  {wavenumber_range[-1]} (cm-1)\n\n")
                input_file.write(f"Npoints {int((wavenumber_range[-1]+1)/spectral_resolution)}\n\n\n")
                input_file.write("absorption\n")
                input_file.write("voigt\n\n")
                input_file.write(f"pressure {pressure:.7f} (bar)\n\n")
                input_file.write("mass 18.\n")
                input_file.write("offset 25.0\n")
                input_file.write("species\n")
                input_file.write(" air   gamma 0.07  n 0.50 t0 296.0 file "+path_data+"1H2-16O__air_a0.broad  delta 0.000\n")
                input_file.write("end\n\n")
                input_file.write("Ncache 1000000\n")
                input_file.write("Nprocs 16\n\n")
                input_file.write(f"Output 1H2-16O_{temperature}K_Voigt_P{pressure}\n\n")
                input_file.write("pffile "+path_data+"1H2-16O__POKAZATEL.pf\n")
                input_file.write("States       "+path_data+"1H2-16O__POKAZATEL.states\n\n")
                input_file.write("Transitions\n")
                for k in range (412): # POKAZATEL goes from 0 to 41200 cm-1
                    lwn = k*100
                    uwn = (k+1)*100
                    if uwn==100:
                        input_file.write("  "+path_data+"1H2-16O__POKAZATEL__00000-00100.trans\n")
                    elif 100<uwn<1000:
                        input_file.write("  "+path_data+"1H2-16O__POKAZATEL__00" +str(lwn)+ "-00"+str(uwn)+".trans\n")
                    elif uwn==1000:
                        input_file.write("  "+path_data+"1H2-16O__POKAZATEL__00900-01000.trans\n")

                    elif 1000<uwn<10000:
                        input_file.write("  "+path_data+"1H2-16O__POKAZATEL__0" +str(lwn)+ "-0"+str(uwn)+".trans\n")

                    elif uwn==10000:
                        input_file.write("  "+path_data+"1H2-16O__POKAZATEL__09900-10000.trans\n")
                    else :
                        input_file.write("  "+path_data+"1H2-16O__POKAZATEL__" +str(lwn)+ "-"+str(uwn)+".trans\n")
                input_file.write("end\n")
            
            try:
                # Run the exocross.exe command with the generated input file
                command = path_exocross+f"xcross.exe < {input_file_name} > {output_file_name}"
                result = subprocess.run(
                    command,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,  # Capture output as text
                    check=True  # Raise an exception if the command returns a non-zero exit code
                )

                # Print the captured standard output (if any)
                if result.stdout:
                    print(result.stdout)

            except subprocess.CalledProcessError as e:
                # Handle errors, including printing the captured standard error
                print(f"Error running xcross.exe: {e}")
                print("Error message:")
                print(e.stderr)
                break
        
            # Optionally, delete the input file after running the command
            if clean:
                os.remove(input_file_name)

            print(f"Processed: {input_file_name}")

    print("All ExoCross runs completed.")