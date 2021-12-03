import os


#os.chdir("/Users/ozgesurer/Desktop/Frescox/experiment/")

def generate_input_file(parameter_values):
    file = 'frescox_inputs/48Ca_template.in'

    with open(file) as f:
        content = f.readlines()
        
    no_p = 0;
    for idx, line in enumerate(content):
        if 'XXXXX' in line:
            no_param = line.count('XXXXX')
            line_temp = line
            for i in range(no_param):
                line_temp = line_temp.replace("XXXXX", str(parameter_values[no_p]), 1) 
                no_p += 1
    
            content[idx] = line_temp
    
    f = open("frescox_inputs/frescox_temp_input.in", "a")
    f.writelines(content)
    f.close()        

def frescox_output(input_file='frescox_inputs/frescox_temp_input.in',
                   output_file='frescox_outputs/48Ca_temp.out'):

    #os.system("../source/frescox < frescox_temp_input.in > 48Ca_temp.out")
    os.system("frescox < frescox_inputs/frescox_temp_input.in > frescox_outputs/48Ca_temp.out")
    with open(output_file) as f:
        content = f.readlines()
        
    sigma_omega_ratio = [] 
    for idline, line in enumerate(content):
        if ('X-S' in line):
            sigma_omega_ratio.append(float(line.split()[4]))
        
    os.remove(input_file)
    os.remove(output_file)
    
    return sigma_omega_ratio

