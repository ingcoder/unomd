# Create trajectory image
import mdtraj as md

def image_molecules(trajectory_input, topology_input, image_output, fraction_of_frames=1.0):
    print("load traj")
    traj = md.load_dcd(trajectory_input, top=topology_input)
    # Slice the trajectory to take only X% of frames. Optional.
    subset = traj[:int(traj.n_frames * fraction_of_frames)]

    print("Start creating image")
    try:
        subset.image_molecules(inplace=True)  # This re-wraps or images the molecules
    except Exception as e:
        print(e)

    try:
        subset.save_dcd(image_output)  # Save the processed trajectory
        print(f"Done! File saved to {image_output}")
    except Exception as e:
        print(e)