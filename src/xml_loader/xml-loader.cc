#include <yafray_config.h>
#include <cstdlib>
#include <cctype>
#include <algorithm>

#ifdef WIN32
	#include <windows.h>
#endif

#include <core_api/scene.h>
#include <core_api/environment.h>
#include <core_api/integrator.h>
#include <core_api/imagefilm.h>
#include <yafraycore/xmlparser.h>
#include <yaf_revision.h>
#include <utilities/console_utils.h>
#include <yafraycore/imageOutput.h>

#include <gui/yafqtapi.h>

using namespace::yafaray;

#include <yafraycore/triangle.h>
#include <yafraycore/meshtypes.h>
#include <limits>

void test(scene_t *scene) {
	int max_to_generate = 320000;
	struct ptd {
		point3d_t p;
		bool t;
	};
	std::vector<ptd> points;
	points.resize(max_to_generate);

	struct disk
	{
		point3d_t p;
		int n_idx;
		float r;
		disk(const point3d_t &p, int n_idx, float r) : p(p), n_idx(n_idx), r(r) {}
	};

	std::vector<vector3d_t> normals;
	std::vector<disk> disks;

	int num_tri_points = 10;
	std::vector<point3d_t> tri_points;
	tri_points.resize(num_tri_points);

	//scene_t::objDataArray_t &meshes = scene->getMeshes();
	scene_t::objDataArray_t &meshes = scene->getMeshes();
	scene_t::objDataArray_t::iterator itr = meshes.begin();
	for(int orig_mesh_count = meshes.size(); orig_mesh_count--; ++itr) {
		triangleObject_t *obj = itr->second.obj;
		int nr_mesh_primitives = obj->numPrimitives();
		const triangle_t **tris = new const triangle_t*[nr_mesh_primitives];
		itr->second.obj->getPrimitives(tris);

		for(int i = 0; i < nr_mesh_primitives; ++i) {
			const triangle_t *t = tris[i];
			point3d_t a, b, c;
			t->getVertices(a, b, c);
			vector3d_t n = t->getNormal();
			normals.push_back(n);
			
			vector3d_t ab = b - a, ac = c - a;
			float area = 0.5 * (ab ^ ac).length();

			float r = 0.1;
			float area_cover = 2*sqrt(2.0);
			int nr_to_keep = std::max((int)(area * area_cover / (r*r*M_PI)), 1);
			int nr_to_generate = nr_to_keep * 10;
			printf("radius: %f\ngenerating: %d/%d\n", r, nr_to_keep, nr_to_generate);
			
			printf("------------\n");
			for(int j = 0; j < nr_to_generate; j++)
			{
				float u = ourRandom();
				float v = ourRandom();
				points[j].p = a + v * ac + (1-v) * u * ab;
				points[j].t = false;
				//printf("(%f,%f,%f)", points[j].p.x, points[j].p.y, points[j].p.z);
			}

			for(int j = 0; j < nr_to_keep; j++)
			{
				int ppoz = 0;
				
				// get the non taken point whose minimum distance 
				// to all taken points is maximum
				float max_mind = std::numeric_limits<float>::min();

				for(int k = 0; k < nr_to_generate; k++)
				{
					float mind = std::numeric_limits<float>::max();
					if(!points[k].t)
					{
						
						// TODO: avoid computing the distances more than once
						for(int l = 0; l < nr_to_generate; l++) {
						//for(int l = 0; l < j; l++) {
						
						//	float d = (points[k].p - disks[l].p).lengthSqr();
							if(!points[l].t) continue;
							float d = (points[k].p - points[l].p).lengthSqr();
							if(d < mind) {
								mind = d;
							}
						}

						if(mind > max_mind) {
							max_mind = mind;
							ppoz = k;
						}
					}
				}

				point3d_t &p = points[ppoz].p;
				points[ppoz].t = true;

				disks.push_back(disk(p, i, r));

				bool hasOrco = false;
				bool hasUV = false;
				int type = 0;
				if(!scene->startGeometry()) exit(0);
				scene->startTriMesh(scene->getNextFreeID(), num_tri_points + 1, num_tri_points, hasOrco, hasUV, type);
				scene->addVertex(p);

				float nz_sq = n.z * n.z;
				int cmap[3] = {0,1,2};
				if(nz_sq < 1e-5) {
					if(n.x * n.x < 1e-5) {
						cmap[1] = 2;
						cmap[2] = 1;
					} else {
						cmap[0] = 2;
						cmap[2] = 0;
					}
					nz_sq = n[cmap[2]] * n[cmap[2]];
				}
				float nx = n[cmap[0]], ny = n[cmap[1]];

				for(int k = 0; k < num_tri_points; k++) {
					float beta_k = k * M_2PI / num_tri_points;
					float cos_beta_k = cos(beta_k);
					float sin_beta_k = sin(beta_k);
					
					float sak_p = (nx * cos_beta_k + ny * sin_beta_k);
					float sinsq_alfa_k = nz_sq / ( sak_p * sak_p + nz_sq);

					float sin_alfa_k = sqrt(sinsq_alfa_k);
					float cos_alfa_k = sqrt(1 - sinsq_alfa_k);
					float r_sin_alfa_k = r * sin_alfa_k;

					tri_points[k][cmap[0]] = p[cmap[0]] + r_sin_alfa_k * cos_beta_k;
					tri_points[k][cmap[1]] = p[cmap[1]] + r_sin_alfa_k * sin_beta_k;
					tri_points[k][cmap[2]] = p[cmap[2]] + r * cos_alfa_k;

					scene->addVertex(tri_points[k]);
				}

				for(int k = 0; k < num_tri_points; k++) {
					scene->addTriangle(0, 1 + k, 1 + (k+1) % num_tri_points, t->getMaterial());
				}
				if(!scene->endTriMesh()) exit(0);
				if(!scene->endGeometry()) exit(0);
			}
		}
		for(std::vector<point3d_t>::iterator p_itr = itr->second.points.begin(); p_itr != itr->second.points.end(); ++p_itr)
			(*p_itr) += point3d_t(3,3,3);
	}
}

int main(int argc, char *argv[])
{
	std::string xmlLoaderVersion = "YafaRay XML loader version 0.2";

	cliParser_t parse(argc, argv, 2, 1, "You need to set at least a yafaray's valid XML file.");

	parse.setAppName(xmlLoaderVersion,
	"[OPTIONS]... <input xml file> [output filename]\n<input xml file> : A valid yafaray XML file\n[output filename] : The filename of the rendered image without extension.\n*Note: If output filename is ommited the name \"yafaray\" will be used instead.");
	
	parse.setOption("pp","plugin-path", false, "Path to load plugins.");
	parse.setOption("vl","verbosity-level", false, "Set verbosity level, options are:\n                                       0 - MUTE (Prints nothing)\n                                       1 - ERROR (Prints only errors)\n                                       2 - WARNING (Prints only errors and warnings)\n                                       3 - INFO (Prints all messages)\n");
	parse.parseCommandLine();
	
#ifdef RELEASE
	std::string version = std::string(VERSION);
#else
	std::string version = std::string(YAF_SVN_REV);
#endif

	renderEnvironment_t *env = new renderEnvironment_t();
	
	// Plugin load
	std::string ppath = parse.getOptionString("pp");
	int verbLevel = parse.getOptionInteger("vl");
	
	if(verbLevel >= 0) yafout.setMasterVerbosity(verbLevel);
	
	if(ppath.empty()) env->getPluginPath(ppath);
	
	if (!ppath.empty())
	{
		Y_INFO << "The plugin path is: " << ppath << yendl;
		env->loadPlugins(ppath);
	}
	else
	{
		Y_ERROR << "Getting plugin path from render environment failed!" << yendl;
		return 1;
	}
	
	std::vector<std::string> formats = env->listImageHandlers();
	
	std::string formatString = "";
	for(size_t i = 0; i < formats.size(); i++)
	{
		formatString.append("                                       " + formats[i]);
		if(i < formats.size() - 1) formatString.append("\n");
	}

	parse.setOption("v","version", true, "Displays this program's version.");
	parse.setOption("h","help", true, "Displays this help text.");
	parse.setOption("op","output-path", false, "Uses the path in <value> as rendered image output path.");
	parse.setOption("f","format", false, "Sets the output image format, available formats are:\n\n" + formatString + "\n                                       Default: tga.\n");
	parse.setOption("t","threads", false, "Overrides threads setting on the XML file, for auto selection use -1.");
	parse.setOption("a","with-alpha", true, "Enables saving the image with alpha channel.");
	parse.setOption("dp","draw-params", true, "Enables saving the image with a settings badge.");
	parse.setOption("ndp","no-draw-params", true, "Disables saving the image with a settings badge (warning: this overrides --draw-params setting).");
	parse.setOption("cs","custom-string", false, "Sets the custom string to be used on the settings badge.");
	parse.setOption("z","z-buffer", true, "Enables the rendering of the depth map (Z-Buffer) (this flag overrides XML setting).");
	parse.setOption("nz","no-z-buffer", true, "Disables the rendering of the depth map (Z-Buffer) (this flag overrides XML setting).");
	
	bool parseOk = parse.parseCommandLine();
	
	if(parse.getFlag("h"))
	{
		parse.printUsage();
		return 0;
	}
	
	if(parse.getFlag("v"))
	{
		Y_INFO << xmlLoaderVersion << yendl << "Built with YafaRay version " << version << yendl;
		return 0;
	}
	
	if(!parseOk)
	{
		parse.printError();
		parse.printUsage();
		return 0;
	}
	
	bool alpha = parse.getFlag("a");
	std::string format = parse.getOptionString("f");
	std::string outputPath = parse.getOptionString("op");
	int threads = parse.getOptionInteger("t");
	bool drawparams = parse.getFlag("dp");
	bool nodrawparams = parse.getFlag("ndp");
	std::string customString = parse.getOptionString("cs");
	bool zbuf = parse.getFlag("z");
	bool nozbuf = parse.getFlag("nz");
	
	if(format.empty()) format = "tga";
	bool formatValid = false;
	
	for(size_t i = 0; i < formats.size(); i++)
	{
		if(formats[i].find(format) != std::string::npos) formatValid = true;
	}
	
	if(!formatValid)
	{
		Y_ERROR << "Couldn't find any valid image format, image handlers missing?" << yendl;
		return 1;
	}
	
	const std::vector<std::string> files = parse.getCleanArgs();
	
	if(files.size() == 0)
	{
		return 0;
	}
	
	std::string outName = "yafray." + format;
	
	if(files.size() > 1) outName = files[1] + "." + format;
	
	std::string xmlFile = files[0];
	
	//env->Debug = debug; //disabled until proper debugging messages are set throughout the core

	// Set the full output path with filename
	if (outputPath.empty())
	{
		outputPath = outName;
	}
	else if (outputPath.at(outputPath.length() - 1) == '/')
	{
		outputPath += outName;
	}
	else if (outputPath.at(outputPath.length() - 1) != '/')
	{
		outputPath += "/" + outName;
	}
	
	scene_t *scene = new scene_t();
	env->setScene(scene);
	paraMap_t render;
	
	bool success = parse_xml_file(xmlFile.c_str(), scene, env, render);
	if(!success) exit(1);

	test(scene);
	
	int width=320, height=240;
	int bx = 0, by = 0;
	render.getParam("width", width); // width of rendered image
	render.getParam("height", height); // height of rendered image
	render.getParam("xstart", bx); // border render x start
	render.getParam("ystart", by); // border render y start
	
	if(threads >= -1) render["threads"] = threads;
	
	if(drawparams)
	{
		render["drawParams"] = true;
		if(!customString.empty()) render["customString"] = customString;
	}
	
	if(nodrawparams) render["drawParams"] = false;
	
	if(zbuf) render["z_channel"] = true;
	if(nozbuf) render["z_channel"] = false;
	
	bool use_zbuf = false;
	render.getParam("z_channel", use_zbuf);
	
	// create output
	colorOutput_t *out = NULL;

	paraMap_t ihParams;
	ihParams["type"] = format;
	ihParams["width"] = width;
	ihParams["height"] = height;
	ihParams["alpha_channel"] = alpha;
	ihParams["z_channel"] = use_zbuf;
	
	imageHandler_t *ih = env->createImageHandler("outFile", ihParams);

	if(ih)
	{
		out = new imageOutput_t(ih, outputPath, bx, by);
		if(!out) return 1;				
	}
	else return 1;
	
	if(! env->setupScene(*scene, render, *out) ) return 1;
	
	scene->render();
	env->clearAll();

	imageFilm_t *film = scene->getImageFilm();

	delete film;
	delete out;
	
	return 0;
}
