# conda run -n omnipose --live-stream python pipeline.py

from functions_pipeline import * 

print("""
Pipeline - SuperSegger-Omnipose

Chose the code to run: 
0 - One FOV Analysis
1 - Multiple FOV's Analysis 
2 - Multiple Experiments Analysis 
3 - To see the documentation
""")

experiment2run = int(input("Pipeline to run: "))

if experiment2run == 3:
    print(help_pipeline.__doc__)
    experiment2run = int(input("Pipeline to run: "))

folder = input("Set the folder to run the analysis: ")
format_image = int(input("Format images (0 = No, 1 = Yes): "))
align = int(input("Align images to the first (0 = No, 1 = Yes): "))
model = int(input("Which model do you want to use? (0 = bact_phase_omni, 1 = bact_fluor_omni): "))

if experiment2run == 0:
    one_fov(folder,format_image,align,model)
elif experiment2run == 1:
    multiple_fovs(folder,format_image,align,model)
elif experiment2run == 2:
    multiple_experiments(folder,format_image,align,model)