import imageio
import subprocess

def make_vid(path = '/Users/alexhill/Desktop/VTK_Ipioni_Comparison/Multi_Panel/Panel_Theta/videos/ppcm_1/', movie_name = 'movie'):

    # Define the path to the directory containing the PNG files
    path_to_dir = path

    # Use subprocess to execute the 'ls' command and capture the output
    #output = subprocess.check_output(["ls", path_to_dir])

    # Define the command to run
    command = ["ls", "-t", path_to_dir]

    # Run the command using subprocess.run
    result = subprocess.run(command, stdout=subprocess.PIPE)

    # Decode the result and print it
    output = result.stdout

    # Split the output into separate file names and store them in a list
    file_names = output.decode().split()

    # Filter the list to include only PNG files
    input_files = [file_name for file_name in file_names if file_name.endswith(".png")]

    #input_files.reverse()

    images = []
    for filename in input_files:
        images.append(imageio.imread(path_to_dir+filename))

    imageio.mimsave(movie_name+'.gif', images[::-1])

    import moviepy.editor as mp

    clip = mp.VideoFileClip(movie_name+".gif")
    clip.write_videofile(movie_name+".mp4")


make_vid('/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs/TSPs/', 'TSPs')

