import numpy as np

def read_coordinates(filename):
    """Reads coordinates from a file."""
    with open(filename, 'r') as f:
        points = []
        for line in f:
            points.append(list(map(float, line.split())))
    return np.array(points)

def write_coordinates(filename, coordinates):
    """Writes coordinates to a file."""
    with open(filename, 'w') as f:
        for point in coordinates:
            f.write(" ".join(f"{x:.6f}" for x in point) + "\n")

def normalize(vector):
    """Normalizes a vector."""
    return vector / np.linalg.norm(vector)

def rodrigues_rotation_matrix(v, z_axis):
    """Computes the rotation matrix using Rodrigues' formula."""
    v = normalize(v)
    k = np.cross(v, z_axis)
    sin_theta = np.linalg.norm(k)
    cos_theta = np.dot(v, z_axis)

    if sin_theta == 0:  # Already aligned
        return np.eye(3)

    k = normalize(k)
    K = np.array([
        [0, -k[2], k[1]],
        [k[2], 0, -k[0]],
        [-k[1], k[0], 0]
    ])
    R = np.eye(3) + sin_theta * K + (1 - cos_theta) * np.dot(K, K)
    return R

def rotate_coordinates(points, R):
    """Applies the rotation matrix to a set of points."""
    return np.dot(points, R.T)

def main():
    input_file = "geom.dat"
    output_file = "geom_rotate.dat"

    # Read points
    points = read_coordinates(input_file)

    # Define vector AD (p4 - p1) and the z-axis
    vector_ad = points[3] - points[0]
    z_axis = np.array([0, 0, 1])

    # Compute rotation matrix
    R = rodrigues_rotation_matrix(vector_ad, z_axis)

    # Rotate all points
    rotated_points = rotate_coordinates(points, R)

    # Write rotated points
    write_coordinates(output_file, rotated_points)

if __name__ == "__main__":
    main()
