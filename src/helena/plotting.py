def plot_geometry(
    ax,
    mesh_coordinates,
    mesh_connections,
    dz_i,
    dr_i,
    Z_mesh_i,
    R_mesh_i,
    image_plotsymmetry,
    image_rotate,
    label_nodes=False,
):
    # Plots mesh geometry onto a 2D plot
    # MeshCoordinates: 2D array of [[radius,height],[radius,height],...]		:: [cm]
    # List of mesh node (vertice) coordinates as floats
    # Expects origin to be top-left corner
    # MeshConnections: 2D array of [[i, j],[i, j],...] 						:: [-]
    # list of indices that identify connections between nodes
    # indices are sequential to MeshCoordinates.

    # Convert strings in scientific notation to floats
    coords = [(float(radius), float(height)) for radius, height in mesh_coordinates]

    # Plot nodes if requested (diagnostic tool)
    if label_nodes:
        for idx, (radius, height) in enumerate(coords, start=1):
            ax.scatter(radius, height, color="red", s=20)
            ax.text(radius, height, f"{idx}", fontsize=8, ha="right", va="bottom")

    # Draw connections
    for n in range(0, len(mesh_connections)):
        # Sequentially read mesh connection indices
        conn = mesh_connections[n]

        # Adjust for zero-based index
        i = int(conn[0]) - 1
        j = int(conn[1]) - 1

        # Extract coordinates at relevant connection indices
        radius1, height1 = coords[i]
        radius2, height2 = coords[j]

        # Determine if line to be plotted is on symmetry axis
        # If so, 'break' and skip to next 'n' to avoid plotting line at symmetry axis
        if image_plotsymmetry:
            # Get max height from current directory [l]
            height_max = round(dz_i * (Z_mesh_i - 1), 2)

            # Check if coordinate one is at the origin
            if radius1 == 0.0 and round(height1, 1) == 0.0:
                # Check if either corresponding height coordinates are at maximum height
                if (
                    abs(height_max - height1) <= dz_i
                    or abs(height_max - height2) <= dz_i
                ):
                    # If both are true, line is on symmetry axis
                    # Skip current line plotting
                    break

            # Check if coordinate two is at the origin
            if radius2 == 0.0 and round(height2, 1) == 0.0:
                # Check if either corresponding height coordinates are at maximum height
                if (
                    abs(height_max - height1) <= dz_i
                    or abs(height_max - height2) <= dz_i
                ):
                    # If both are true, line is on symmetry axis
                    # Skip current line plotting
                    break

        # Only works if rotating 90 to the right
        if image_rotate and image_plotsymmetry:
            radius1, height1 = height1, radius1
            radius2, height2 = height2, radius2
            ax.plot([radius1, radius2], [height1, height2], color="dimgrey", lw=2)
            ax.plot([radius2, radius1], [-height1, -height2], color="dimgrey", lw=2)

        # Only works if rotating 90 to the right
        elif image_rotate and not image_plotsymmetry:
            radius1, height1 = height1, radius1
            radius2, height2 = height2, radius2
            height1, height2 = -height2 + R_mesh_i * dr_i, -height1 + R_mesh_i * dr_i
            ax.plot([radius1, radius2], [height1, height2], color="dimgrey", lw=2)

        # Works for all cases
        elif not image_rotate and image_plotsymmetry:
            ax.plot([radius1, radius2], [height1, height2], color="dimgrey", lw=2)
            ax.plot([-radius2, -radius1], [height1, height2], color="dimgrey", lw=2)

        # if no rotation or symmetry
        else:
            ax.plot([radius1, radius2], [height1, height2], color="dimgrey", lw=2)
