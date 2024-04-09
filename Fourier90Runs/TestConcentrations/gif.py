import imageio.v2 as imageio
import subprocess

def make_vid(path = '/Users/alexhill/Desktop/VTK_Ipioni_Comparison/Multi_Panel/Panel_Theta/videos/ppcm_1/', movie_name = 'movie'):

    # Define the path to the directory containing the PNG files
    path_to_dir = path

    # Use subprocess to execute the 'ls' command and capture the output
    #ßoutput = subprocess.check_output(["ls", path_to_dir])

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
    gif_file = movie_name+'.gif'
    # Output MP4 file name
    mp4_file = movie_name+".mp4"

    # Read the GIF file
    gif = imageio.mimread(gif_file)

    # Write the MP4 file
    with imageio.get_writer(mp4_file, fps=10, macro_block_size = 1) as writer:  # Adjust fps as needed
        for frame in gif:
            writer.append_data(frame)



make_vid(path='/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs/TestConcentrations/sig/', movie_name='sig')
make_vid(path='/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs/TestConcentrations/noise/', movie_name='noise')



make_vid(path = '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs/TestConcentrations/sig/', movie_name = 'sig')
make_vid(path = '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs/TestConcentrations/noise/', movie_name = 'noise')

