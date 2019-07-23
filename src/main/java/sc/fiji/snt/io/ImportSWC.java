package sc.fiji.snt.io;

import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultUndirectedGraph;
import org.scijava.util.ColorRGB;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ImportSWC {
    	/**
	 * Imports SWC data with advanced settings. The SWC format is described in
	 * <a href="https://www.ncbi.nlm.nih.gov/pubmed/9821633">PMID 9821633</a> and
	 * e.g.,
	 * <a href="http://www.neuromorpho.org/myfaq.jsp#qr3">neuromorpho.org</a> It
	 * is named after the initials of Stockley, Wheal, and Cole, who earlier
	 * developed a pioneer system for morphometric reconstructions
	 * (<a href="https://www.ncbi.nlm.nih.gov/pubmed/8321013">PMID 8321013</a>).
	 * <p>
	 * Annoyingly, While the SWC specification details the usage of of world
	 * coordinates in microns, some published SWC files have adopted image (pixel)
	 * coordinates, which is inappropriate and less useful (an example of the
	 * latter seems to part of the DIADEM Challenge data set). In addition, it's
	 * not clear what the "radius" column is meant to mean in such files.
	 * </p>
	 *
	 * @param br the character stream containing the data
	 * @param assumeCoordinatesInVoxels If true, the SWC coordinates are assumed
	 *          to be in image coordinates ("pixels"). Note that in this case,
	 *          radii will be scaled by the minimum voxel separation. This
	 *          workaround seems to be required to properly import unscaled files
	 * @param xOffset the offset to be applied to all X coordinates. May be useful
	 *          to import data obtained from multiple "un-stitched" fields of
	 *          view. Default is 0.
	 * @param yOffset the offset to be applied to all Y coordinates. May be useful
	 *          to import data obtained from multiple "un-stitched" fields of
	 *          view. Default is 0.
	 * @param zOffset the offset to be applied to all Z coordinates. May be useful
	 *          to import data obtained from multiple "un-stitched" fields of
	 *          view. Default is 0.
	 * @param xScale the scaling factor for all X coordinates. Useful to import
	 *          data onto downsampled images. Default is 1.
	 * @param yScale the scaling factor for all Y coordinates. Useful to import
	 *          data onto downsampled images. Default is 1.
	 * @param zScale the scaling factor for all Z coordinates. Useful to import
	 *          data onto downsampled images. Default is 1.
	 * @param replaceAllPaths If true, all existing Paths will be deleted before
	 *          the import. Default is false.
	 * @return true, if import was successful
	 */
	protected Graph<double[], double[]> importSWC(final BufferedReader br,
							  final boolean assumeCoordinatesInVoxels, final double xOffset,
							  final double yOffset, final double zOffset, final double xScale,
							  final double yScale, final double zScale, final boolean replaceAllPaths)
	{
		final Pattern pEmpty = Pattern.compile("^\\s*$");
		final Pattern pComment = Pattern.compile("^([^#]*)#.*$");

		final TreeSet<double[]> nodes = new TreeSet<>((doubles, t1) -> {
			if( doubles[0] > t1[0] ) return 1;
			else if( doubles[0] == t1[0] ) return 0;
			else return -1;
		});

		String line;
		try {
			while ((line = br.readLine()) != null) {
				final Matcher mComment = pComment.matcher(line);
				line = mComment.replaceAll("$1").trim();
				final Matcher mEmpty = pEmpty.matcher(line);
				if (mEmpty.matches()) continue;
				final String[] fields = line.split("\\s+");
				if (fields.length < 7) {
					System.err.println("Wrong number of fields (" + fields.length + ") in line: " + line);
					return null;
				}
				try {
					final int id = Integer.parseInt(fields[0]);
					final int type = Integer.parseInt(fields[1]);
					final double x = xScale * Double.parseDouble(fields[2]) + xOffset;
					final double y = yScale * Double.parseDouble(fields[3]) + yOffset;
					final double z = zScale * Double.parseDouble(fields[4]) + zOffset;
					double radius;
					try {
						radius = Double.parseDouble(fields[5]);
					} catch (final NumberFormatException ignored) {
						radius = 0; // files in which radius is set to NaN
					}
					final int previous = Integer.parseInt(fields[6]);
					nodes.add(new double[]{id, type, x, y, z, radius, previous});
				}
				catch (final NumberFormatException nfe) {
					System.err.println("There was a malformed number in line: " + line);
					return null;
				}
			}
		}
		catch (final IOException exc) {
			System.err.println("IO ERROR"+ exc);
			return null;
		}
		return importNodes(null, nodes, null, assumeCoordinatesInVoxels);
	}

	private Graph<double[], double[]> importNodes(final String descriptor,
	                            final TreeSet<double[]> points, final ColorRGB color,
	                            final boolean assumeCoordinatesInVoxels)
	{
		Graph<double[], double[]> graph = new DefaultUndirectedGraph<double[], double[]>();

		final Map<Integer, double[]> idToSWCPoint = new HashMap<>();
		final List<double[]> primaryPoints = new ArrayList<>();

		final int firstImportedPathIdx = -1;//size();
		for (final double[] point : points) {
			idToSWCPoint.put((int)point[0], point);
			if ((int)point[6] == -1) {
				primaryPoints.add(point);
			}
			else {
				final double[] previousPoint = idToSWCPoint.get((int)point[6]);
				if (previousPoint != null) {
					point[6] = previousPoint[0];
					graph.addEdge(point, previousPoint);
				}
			}
		}

		// TODO reactivate this code
		/*if (spacingIsUnset) {
			if (!headless) SNTUtils.log(
				"Inferring pixel spacing from imported points... ");
			boundingBox.inferSpacing(points);
			x_spacing = boundingBox.xSpacing;
			y_spacing = boundingBox.ySpacing;
			z_spacing = boundingBox.zSpacing;
			// NB: we must leave boundingBox origin unset so that it's dimensions can be properly computed
			spacingIsUnset = false;
		}*/

		double x_spacing = 1;
		double y_spacing = 1;
		double z_spacing = 1;

		// We'll now iterate (again!) through the points to fix some ill-assembled
		// files that do exist in the wild defined in pixel coordinates!
		if (assumeCoordinatesInVoxels) {
			final double minimumVoxelSpacing = Math.min(Math.abs(x_spacing),
					Math.min(Math.abs(y_spacing), Math.abs(z_spacing)));

			final Iterator<double[]> it = points.iterator();
			while (it.hasNext()) {
				final double[] point = it.next();

				point[2] *= x_spacing;
				point[3] *= y_spacing;
				point[4] *= z_spacing;
				// this just seems to be the convention in the broken files we've came
				// across
				point[5] *= minimumVoxelSpacing;

				// If the radius is set to near zero, then artificially set it to half
				// of the voxel spacing so that something* appears in the 3D Viewer!
				if (Math.abs(point[5]) < 0.0000001)
					point[5] = minimumVoxelSpacing / 2;
			}
		}

		// TODO resume here

		// FIXME: This is really slow with large SWC files
		final HashMap<SWCPoint, Path> pointToPath = new HashMap<>();
		final PriorityQueue<SWCPoint> backtrackTo = new PriorityQueue<>(
			primaryPoints);
		final HashMap<Path, SWCPoint> pathStartsOnSWCPoint = new HashMap<>();
		final HashMap<Path, PointInImage> pathStartsAtPointInImage =
			new HashMap<>();

		SWCPoint start;
		Path currentPath;
		while ((start = backtrackTo.poll()) != null) {
			currentPath = new Path(x_spacing, y_spacing, z_spacing, spacing_units);
			currentPath.createCircles();
			if (start.getPreviousPoint() != null) {
				final SWCPoint beforeStart = start.getPreviousPoint();
				pathStartsOnSWCPoint.put(currentPath, beforeStart);
				pathStartsAtPointInImage.put(currentPath, beforeStart);
				currentPath.addNode(beforeStart);
			}

			// Now we can start adding points to the path:
			SWCPoint currentPoint = start;
			while (currentPoint != null) {
				currentPath.addNode(currentPoint);
				pointToPath.put(currentPoint, currentPath);

				if (currentPoint.getNextPoints().size() > 0) {
					final SWCPoint newCurrentPoint = currentPoint.getNextPoints().get(0);
					currentPoint.getNextPoints().remove(0);
					backtrackTo.addAll(currentPoint.getNextPoints());
					currentPoint = newCurrentPoint;
				}
				else {
					currentPath.setSWCType(currentPoint.type); // Assign point
					// type to path
					currentPoint = null;
				}
			}
			currentPath.setGuessedTangents(2);
			addPath(currentPath);
		}

		// Set the start joins:
		for (int i = firstImportedPathIdx; i < size(); i++) {
			final Path p = getPath(i);
			final SWCPoint swcPoint = pathStartsOnSWCPoint.get(p);
			if (descriptor != null) {
				p.setName(p.getName() + "{" + descriptor + "}");
				p.setColorRGB(color);
			}
			if (swcPoint == null) {
				p.setIsPrimary(true);
				continue;
			}
			final Path previousPath = pointToPath.get(swcPoint);
			final PointInImage pointInImage = pathStartsAtPointInImage.get(p);
			p.setStartJoin(previousPath, pointInImage);
		}

		resetListeners(null, true);

		// Infer fields for when an image has not been specified. We'll assume
		// the image dimensions to be those of the coordinates bounding box.
		// This allows us to open a SWC file without a source image
		{
			if (boundingBox == null)
				boundingBox = new BoundingBox();
			boundingBox.compute(((TreeSet<? extends SNTPoint>) points).iterator());

			// If a plugin exists, warn user if its image cannot render imported nodes
			if (plugin != null && plugin.getImagePlus() != null) {

				final BoundingBox pluginBoundingBox = new BoundingBox();
				pluginBoundingBox.setOrigin(new PointInImage(0, 0, 0));
				pluginBoundingBox.setDimensions(plugin.width, plugin.height, plugin.depth);
				if (!pluginBoundingBox.contains(boundingBox)) {
					SNTUtils.warn("Some nodes lay outside the image volume: you may need to "
							+ "adjust import options or resize current image canvas");
				}
			}
		}
		return true;
	}
}
