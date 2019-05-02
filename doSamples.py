import os
import subprocess

#'biological'
for root, dir, files in os.walk('reduced_samples_V_100to1000_k_0.5to10V'):
    for item in files:
        inputpath = os.path.join(root, item)
        outputpath = 'solved_' + os.path.splitext(inputpath)[0]
        reducedout = 'reduced_' + os.path.splitext(inputpath)[0] + '.cm'
        outputdir = os.path.split(outputpath)[0]
        reduceddir = os.path.split(reducedout)[0]
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        if not os.path.exists(reduceddir):
            os.makedirs(reduceddir)
        try:
            haveToDo = True
            try:
                with open(outputpath + '.out', 'r') as f:
                    if f.read() != "TIMEOUT":
                        haveToDo = False
                        print("Skipped {}, because it has already been solved".format(outputpath))
                    else:
                        print("{} was previously timed out, retrying...".format(outputpath))
            except FileNotFoundError as e:
                haveToDO = True
            if haveToDo:
                result = subprocess.run(["yoshiko/build/yoshiko.2.2.0", "-f", inputpath, "-v", "2", "-o", outputpath, "-threads", "1", "-T", "300"], stdout=subprocess.PIPE, timeout=500) # "-H", "1"
                # result = subprocess.run(["yoshiko/build/yoshiko.2.2.0", "-f", inputpath, "-v", "2", "-o", outputpath, "-threads", "1", "-save-reduced", reducedout], stdout=subprocess.PIPE, timeout=500)
                with open(outputpath + '.out', 'wb') as f:
                    f.write(result.stdout)
        except subprocess.TimeoutExpired as e:
            print("Timeout on {}".format(outputpath))
            with open(outputpath + '.out', 'w') as f:
                    f.write("TIMEOUT")