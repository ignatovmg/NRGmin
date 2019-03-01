#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <jansson.h>
#include <getopt.h>

#include "mol2/json.h"
#include "mol2/benergy.h"
#include "mol2/gbsa.h"
#include "mol2/icharmm.h"
#include "mol2/minimize.h"
#include "mol2/nbenergy.h"
#include "mol2/pdb.h"

#define __TOL__ 5E-4

#define ERR_MSG(fmt, ...) {                                   \
    fprintf(stderr, fmt "in file %s, function %s, line %i\n"  \
            "Exiting ...\n", ##__VA_ARGS__, __FILE__,         \
            __func__, __LINE__);                              \
    exit(EXIT_FAILURE);                                       \
}

#define MYFOPEN(fp, path, spec)  fopen(path, spec); {         \
    if (fp == NULL) {                                         \
        ERR_MSG("ERROR opening file %s\n", path)              \
    }                                                         \
}

#define MYCALLOC(ptr, n, size)  calloc(n, size); {            \
    if (ptr == NULL) {                                        \
        ERR_MSG("ERROR allocating %i bytes for variable %s\n",\
            (int)((n))*(int)((size)), #ptr)                   \
    }                                                         \
}

#define MAGIC_ARGS(arg) {                                     \
	char name[] = #arg;                                       \
	char* pos = strchr(name, '_');                            \
	if (pos != NULL) *pos = '-';                              \
	if (strcmp(name, long_options[option_index].name) == 0) { \
			arg = optarg;                                     \
	}                                                         \
}

struct energy_prm {

	struct mol_atom_group *ag;
	struct agsetup  *ag_setup;
	struct acesetup *ace_setup;
    struct springset_pairs* sprst_pairs;
    struct springset_points* sprst_points;
};


struct pair_spring
{
        struct mol_atom_group *ag;  /**< affected atomgroup */
        int    laspr[2];      /**< list of atoms */
        double lnspr;
        double erspr;
        double fkspr;      /**< force constant */
};

struct point_spring
{
        struct mol_atom_group *ag;  /**< affected atomgroup */
        int    naspr;      /**< number of affected atoms */
        int   *laspr;      /**< list of atoms */
        //double lnspr;
        double fkspr;      /**< force constant */
        double X0, Y0, Z0; /**< anchor point */
};

struct springset_pairs
{
        int nsprings;       /**< number of springs */
        struct pair_spring *springs;  /**< array of springs */
};

struct springset_points
{
        int nsprings;       /**< number of springs */
        struct point_spring *springs;  /**< array of springs */
};


void read_fix(char *ffile, int *nfix, size_t **fix);

struct springset_pairs* read_springset_pairs(struct mol_atom_group* ag, char *sfile);

struct springset_points* read_springset_points(struct mol_atom_group* ag, char *sfile);

void free_springset_pairs(struct springset_pairs** sprst);

void free_springset_points(struct springset_points** sprst);

void springeng_pair(struct springset_pairs *sprst, double* een);

void springeng_point(struct springset_points *sprst, double* een);

static lbfgsfloatval_t energy_func(
	void* restrict prm,
	const double* restrict array,
	double* restrict gradient,
	const int array_size,
	const lbfgsfloatval_t step);
	
static void fprint_energy_terms(FILE* stream, void* restrict prm, char* prefix);
	
void usage_message(char** argv);

void help_message(void);

int main(int argc, char** argv) 
{
	//static int verbose_flag = 0;
	static int ace_flag = 0;
	char* protocol = NULL;
	char* pdb     = NULL;
	char* psf     = NULL;
	char* prm     = NULL;
	char* rtf     = NULL;
	char* out     = NULL;
	char* rec_pdb = NULL;
	char* rec_psf = NULL;
	char* lig_pdb = NULL;
	char* lig_psf = NULL;
	char* lig_jsn = NULL;
	int   nsteps  = 1000;
	
	static struct option long_options[] =
	{
		//{"verbose", no_argument, &verbose_flag, 1},
		{"rtf",  required_argument, 0, 0},
		{"prm",  required_argument, 0, 0},
		{"out",  required_argument, 0, 0},
		{"protocol", required_argument, 0, 0},
		{"pdb", required_argument, 0, 0},
		{"psf", required_argument, 0, 0},
		{"rec-pdb",  required_argument, 0, 0},
		{"rec-psf",  required_argument, 0, 0},
		{"lig-pdb",  required_argument, 0, 0},
		{"lig-psf",  required_argument, 0, 0},
		{"lig-jsn",  required_argument, 0, 0},
		{"nsteps",  required_argument, 0, 0},
		{"gbsa", no_argument, &ace_flag, 1},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	};

	int option_index = 0;
	int opt;
	while (1) {
		opt = getopt_long_only (argc, argv, "h", long_options, &option_index);
		
		
		if (opt == -1) {
			break;
		}
		
		switch (opt) {	
			case 0:
				if (long_options[option_index].flag != 0)
					break;
				printf ("Option %s", long_options[option_index].name);
				if (optarg)
					printf (" = %s", optarg);
				printf ("\n");
				break;
			
			case 'h':
				usage_message(argv);
				exit(EXIT_SUCCESS);
				break;
				
			case '?':
				usage_message(argv);
				exit(EXIT_FAILURE);
				break;
				
			default:
				usage_message(argv);
				exit(EXIT_FAILURE);
				break;
		}
		
		MAGIC_ARGS(rtf);
		MAGIC_ARGS(prm);
		MAGIC_ARGS(out);
		MAGIC_ARGS(pdb);
		MAGIC_ARGS(psf);
		MAGIC_ARGS(protocol);
		MAGIC_ARGS(rec_pdb);
		MAGIC_ARGS(rec_psf);
		MAGIC_ARGS(lig_pdb);
		MAGIC_ARGS(lig_psf);
		MAGIC_ARGS(lig_jsn);
		
		if (strcmp("nsteps", long_options[option_index].name) == 0) {
			nsteps = atoi(optarg);
		}
		
	}
	
	struct mol_atom_group *ag = NULL;
	struct mol_atom_group_list* aglist = NULL;
	
	if (rtf == NULL || prm == NULL) {
		fprintf(stderr, "RTF and PRM are required\n");
		usage_message(argv);
		exit(EXIT_FAILURE);
	}
	
	if (out == NULL) {
		out = "min.pdb";
	}
	
	
	if (pdb != NULL || psf != NULL ) {
		if (pdb == NULL) {
			fprintf(stderr, "Both PDB and PSF must be specified\n");
			usage_message(argv);
			exit(EXIT_FAILURE);
		}
		if (psf == NULL) {
			fprintf(stderr, "Both PDB and PSF must be specified\n");
			usage_message(argv);
			exit(EXIT_FAILURE);
		}
		if (rec_pdb != NULL || rec_psf != NULL) {
			fprintf(stderr, "Full molecule PDB and PSF are already specified. Receptor PDB and PSF are redundant\n");
			usage_message(argv);
			exit(EXIT_FAILURE);
		}
		if (lig_pdb != NULL || lig_psf != NULL || lig_jsn != NULL) {
			fprintf(stderr, "Full molecule PDB and PSF are already specified. Ligand PDB and PSF or JSON are redundant\n");
			usage_message(argv);
			exit(EXIT_FAILURE);
		}
		
		//ag = mol_read_pdb(pdb);
		aglist = mol_read_pdb_models(pdb);
		for (size_t i = 0; i < aglist->size; i++) {
			mol_atom_group_read_geometry(&aglist->members[i], psf, prm, rtf);
		}
		
	} else {
		struct mol_atom_group *rec_ag, *lig_ag;
	
		if (rec_pdb == NULL || rec_psf == NULL) {
			fprintf(stderr, "Please provide receptor PDB and PSF\n");
			usage_message(argv);
			exit(EXIT_FAILURE);
		} else {
		
			rec_ag = mol_read_pdb(rec_pdb);
			mol_atom_group_read_geometry(rec_ag, rec_psf, prm, rtf);
		}
	
		if (lig_pdb != NULL || lig_psf != NULL) {
			if (lig_pdb == NULL || lig_psf == NULL) {
				fprintf(stderr, "Both ligand pdb and psf must be specified\n");
				usage_message(argv);
				exit(EXIT_FAILURE);
			}
			if (lig_jsn != NULL) {
				fprintf(stderr, "Ligand PDB and PSF are already specified. Ligand JSON is redundant\n");
				usage_message(argv);
				exit(EXIT_FAILURE);
			}
		
			lig_ag = mol_read_pdb(lig_pdb);
			mol_atom_group_read_geometry(lig_ag, lig_psf, prm, rtf);
			
		} else {
			if (lig_jsn == NULL) {
				fprintf(stderr, "Please provide ligand PDB and PSF or JSON\n");
				usage_message(argv);
				exit(EXIT_FAILURE);
			}
		
			lig_ag = mol_read_json(lig_jsn);
		}
		
		ag = mol_atom_group_join(rec_ag, lig_ag);
		mol_atom_group_free(lig_ag);
		mol_atom_group_free(rec_ag);
	}

	FILE* outfile = MYFOPEN(outfile, out, "w");

	for (int modeli = 0; modeli < aglist->size; modeli++) {
		fprintf(outfile, "MODEL %i\n", (modeli + 1));

		ag = &aglist->members[modeli];
		ag->gradients = MYCALLOC(ag->gradients, ag->natoms, sizeof(struct mol_vector3));
		mol_fixed_init(ag);	
		mol_fixed_update(ag, 0, NULL);
		
		struct agsetup ags;
		init_nblst(ag, &ags);
		update_nblst(ag, &ags);
		
		struct acesetup ace_setup;
		if (ace_flag == 1) {
			ace_setup.efac = 0.5;
			ace_ini(ag, &ace_setup);
			ace_fixedupdate(ag, &ags, &ace_setup);
			ace_updatenblst(&ags, &ace_setup);
		}
		
		struct springset_pairs *sprst_pairs = NULL;
		struct springset_points *sprst_points = NULL;      
		
		/*if (protocol != NULL && nsteps != -1) {
			fprintf(stderr, "Please provide either protocol XOR nsteps\n");
			exit(EXIT_FAILURE);
		}*/
		
		struct energy_prm engpar;
		       
		if (protocol != NULL) {
			FILE* prot_file = MYFOPEN(prot_file, protocol, "r");
			int  cur_nsteps;
			char fix_path[1024];
			char spr_pair_path[1024];
			char spr_point_path[1024];
			
			int c;
			while ((c = fscanf(prot_file, "%i %s %s %s", &cur_nsteps, fix_path, 
					      spr_pair_path, spr_point_path)) != EOF) {
				if (c != 4) { 
					ERR_MSG("Wrong protocol file format (%i words read)\n", c); 
				}
				
				if (cur_nsteps < 0) { 
					ERR_MSG("Wrong protocol file format (number of steps must be non-negative)\n"); 
				}

				int     nfix = 0;
				size_t* fix  = NULL;

				if (strcmp(fix_path, ".") != 0) {
					read_fix(fix_path, &nfix, &fix);
					mol_fixed_update(ag, nfix, fix);
					update_nblst(ag, &ags);

					if (ace_flag == 1) {
						ace_fixedupdate(ag, &ags, &ace_setup);
						ace_updatenblst(&ags, &ace_setup);
					}
						
					free(fix);
					fix = NULL;
				} 

				if (strcmp(spr_pair_path, ".") != 0) {
					sprst_pairs = read_springset_pairs(ag, spr_pair_path);
				}

				if (strcmp(spr_point_path, ".") != 0) {
					sprst_points = read_springset_points(ag, spr_point_path);
				}
				    
				engpar.ag = ag;
				engpar.ag_setup = &ags;
				engpar.sprst_pairs = sprst_pairs;
				engpar.sprst_points = sprst_points;

				if (ace_flag == 1) {
					engpar.ace_setup = &ace_setup;
				} else {
					engpar.ace_setup = NULL;
				}

				//fprintf(outfile, "REMARK Starting energies:\n");
				fprint_energy_terms(outfile, &engpar, "REMARK START ");

				if (cur_nsteps > 0) {
					mol_minimize_ag(MOL_LBFGS, cur_nsteps, __TOL__, ag, (void *)(&engpar), energy_func);
				}
			}
		} else {
			if (nsteps < 0) {
				ERR_MSG("Number of steps must be non-negative (nsteps = %i)\n", nsteps); 
			}
			
			engpar.ag = ag;
			engpar.ag_setup  = &ags;
			engpar.sprst_pairs = NULL;
			engpar.sprst_points = NULL;
			if (ace_flag == 1) {
				engpar.ace_setup = &ace_setup;
			} else {
				engpar.ace_setup = NULL;
			}
			
			//fprintf(outfile, "REMARK Starting energies:\n");
			fprint_energy_terms(outfile, &engpar, "REMARK START ");
			
			if (nsteps > 0) {
				mol_minimize_ag(MOL_LBFGS, nsteps, __TOL__, ag, (void *)(&engpar), energy_func);
			}
			
		}

		//printf("\nFinal energies:\n");
		fprintf(outfile, "REMARK\n");		
		fprint_energy_terms(outfile, &engpar, "REMARK FINAL ");
		mol_fwrite_pdb(outfile, ag);
		fprintf(outfile, "ENDMDL\n");
		fflush(outfile);	

		free_springset_pairs(&sprst_pairs);
		free_springset_points(&sprst_points);
	}

	fclose(outfile);	
	mol_atom_group_list_free(aglist);
	
	return EXIT_SUCCESS;
}


static lbfgsfloatval_t energy_func(
	void* restrict prm,
	const double* restrict array,
	double* restrict gradient,
	const int array_size,
	const lbfgsfloatval_t step)
{
	lbfgsfloatval_t energy = 0.0;
	struct energy_prm* energy_prm = (struct energy_prm*) prm;
	
	if (array != NULL) {
		mol_atom_group_set_actives(energy_prm->ag, array);
	}
	bool updated = check_clusterupdate(energy_prm->ag, energy_prm->ag_setup);
	if (updated) {
		if (energy_prm->ace_setup != NULL) {
			ace_updatenblst(energy_prm->ag_setup, energy_prm->ace_setup);
		}
	}
	
	mol_zero_gradients(energy_prm->ag);
	
	if (energy_prm->ace_setup != NULL) {
		aceeng(energy_prm->ag, &energy, energy_prm->ace_setup, energy_prm->ag_setup);
	}
	
	vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
	vdwengs03(1.0, energy_prm->ag_setup->nblst->nbcof, energy_prm->ag, &energy,
	          energy_prm->ag_setup->nf03, energy_prm->ag_setup->listf03);
	beng(energy_prm->ag, &energy);
	aeng(energy_prm->ag, &energy);
	teng(energy_prm->ag, &energy);
	ieng(energy_prm->ag, &energy);
	
	if (energy_prm->sprst_pairs != NULL) {
		springeng_pair(energy_prm->sprst_pairs, &energy);
	}
	if (energy_prm->sprst_points != NULL) {
		springeng_point(energy_prm->sprst_points, &energy);
	}

	if (gradient != NULL) {
		for (int i = 0; i < array_size / 3; i++ ) {
			int atom_i = energy_prm->ag->active_atoms->members[i];
			gradient[3*i  ] = -energy_prm->ag->gradients[atom_i].X;
			gradient[3*i+1] = -energy_prm->ag->gradients[atom_i].Y;
			gradient[3*i+2] = -energy_prm->ag->gradients[atom_i].Z;
		}
	}
	
	return energy;
}

static void fprint_energy_terms(FILE* stream, void* restrict prm, char* prefix)
	//const double* restrict array,
	//const int array_size,
	//const lbfgsfloatval_t step)
{
	lbfgsfloatval_t energy = 0.0;
	lbfgsfloatval_t total = 0.0;
	struct energy_prm* energy_prm = (struct energy_prm*) prm;
	
	//if (array != NULL) {
	//	mol_atom_group_set_actives(energy_prm->ag, array);
	//}
	bool updated = check_clusterupdate(energy_prm->ag, energy_prm->ag_setup);
	if (updated) {
		if (energy_prm->ace_setup != NULL) {
			ace_updatenblst(energy_prm->ag_setup, energy_prm->ace_setup);
		}
	}

	char fmt[1024];
	strcpy(fmt, prefix);

	if (energy_prm->ace_setup != NULL) {
		aceeng(energy_prm->ag, &energy, energy_prm->ace_setup, energy_prm->ag_setup);
		strcpy(fmt, prefix);
		fprintf(stream, strcat(fmt, "ACE: % .3f\n"), energy);
		total += energy;
		energy = 0.0;
	}
	
	vdweng(energy_prm->ag, &energy, energy_prm->ag_setup->nblst);
	strcpy(fmt, prefix);
	fprintf(stream, strcat(fmt, "VWD: % .3f\n"), energy);
	total += energy;
	energy = 0.0;
	
	vdwengs03(1.0, energy_prm->ag_setup->nblst->nbcof, energy_prm->ag, &energy,
	          energy_prm->ag_setup->nf03, energy_prm->ag_setup->listf03);
	strcpy(fmt, prefix);
	fprintf(stream, strcat(fmt, "VWD03: % .3f\n"), energy);
	total += energy;
	energy = 0.0;
	
	beng(energy_prm->ag, &energy);
	strcpy(fmt, prefix);
	fprintf(stream, strcat(fmt, "Bonded: % .3f\n"), energy);
	total += energy;
	energy = 0.0;
	
	aeng(energy_prm->ag, &energy);
	strcpy(fmt, prefix);
	fprintf(stream, strcat(fmt, "Angles: % .3f\n"), energy);
	total += energy;
	energy = 0.0;
	
	teng(energy_prm->ag, &energy);
	strcpy(fmt, prefix);
	fprintf(stream, strcat(fmt, "Torsions: % .3f\n"), energy);
	total += energy;
	energy = 0.0;
	
	ieng(energy_prm->ag, &energy);
	strcpy(fmt, prefix);
	fprintf(stream, strcat(fmt, "Impropers: % .3f\n"), energy);
	total += energy;
	energy = 0.0;
	
	if (energy_prm->sprst_pairs != NULL) {
		springeng_pair(energy_prm->sprst_pairs, &energy);
		strcpy(fmt, prefix);
		fprintf(stream, strcat(fmt, "Pairsprings: % .3f\n"), energy);
		total += energy;
		energy = 0.0;
	}
	if (energy_prm->sprst_points != NULL) {
		springeng_point(energy_prm->sprst_points, &energy);
		strcpy(fmt, prefix);
		fprintf(stream, strcat(fmt, "Pointsprings: % .3f\n"), energy);
		total += energy;
		energy = 0.0;
	}

	strcpy(fmt, prefix);	
	fprintf(stream, strcat(fmt, "Total: % .3f\n"), total);
}

void read_fix(char *ffile, int *nfix, size_t **fix) 
{ 
	int linesz = 91; 
	char* buffer = (char*)MYCALLOC(buffer, linesz, sizeof(char)); 

	*nfix=0; 
	FILE* fp = MYFOPEN(fp, ffile, "r");
	
	while (fgets(buffer, linesz-1, fp) != NULL) {
		if(!strncmp(buffer,"ATOM",4))(*nfix)++; 
	}

	rewind(fp); 
	*fix = MYCALLOC(*fix, *nfix, sizeof(size_t)); 
	//fp = MYFOPEN(fp, ffile, "r"); 
	int na = 0; 

	while(fgets(buffer, linesz-1, fp) != NULL) { 
		if(!strncmp(buffer,"ATOM",4)) { 
			(*fix)[na] = atoi(buffer+4)-1; 
			na++; 
		} 
	} 

	free(buffer); 
	fclose(fp); 
}

struct springset_pairs* read_springset_pairs(struct mol_atom_group* ag, char *sfile) {
	FILE* fp = MYFOPEN(fp, sfile, "r");
	
	struct springset_pairs* sprst;
	sprst = MYCALLOC(sprst, 1, sizeof(struct springset_pairs));
	
	int c;
	if (fscanf(fp, "%i", &sprst->nsprings) != 1) {
		ERR_MSG("Wrong spring file format\n");
	}
	
	sprst->springs = MYCALLOC(sprst->springs, sprst->nsprings, sizeof(struct pair_spring));
	struct pair_spring* sprs = sprst->springs;
	
	int id = 0;
	char name1[8], name2[8];
	int   aid1,   aid2;
	while (id < sprst->nsprings) {
		//             id1 n1 id2 n2 len err fk
	    c = fscanf(fp, "%lf %lf %lf %i %s %i %s", &sprs[id].lnspr, &sprs[id].erspr, 
	                    &sprs[id].fkspr, &aid1, name1, &aid2, name2);
	               
	    //printf("%lf %lf %lf %i %s %i %s\n", sprs[id].lnspr, sprs[id].erspr, sprs[id].fkspr, aid1, name1, aid2, name2);
		if (c != 7) {
			ERR_MSG("Wrong spring file format\n");
		}
		
		sprs[id].ag = ag;
		//printf("%s %s\n", ag->atom_name[aid1-1], ag->atom_name[aid2-1]);
		if (strstr(ag->atom_name[aid1-1], name1) == NULL || strstr(ag->atom_name[aid2-1], name2) == NULL) {
			ERR_MSG("Inconsistent numbering in file %s\n", sfile);
		}
		
		sprs[id].laspr[0] = aid1 - 1;
		sprs[id].laspr[1] = aid2 - 1;
		
		id++;
	}
	
	fclose(fp);
	return sprst;
}

struct springset_points* read_springset_points(struct mol_atom_group* ag, char *sfile) {
	FILE* fp = MYFOPEN(fp, sfile, "r");
	
	struct springset_points* sprst;
	sprst = MYCALLOC(sprst, 1, sizeof(struct springset_points));
	
	int c;
	if (fscanf(fp, "%i", &sprst->nsprings) != 1) {
		ERR_MSG("Wrong spring file format\n");
	}
	
	sprst->springs = MYCALLOC(sprst->springs, sprst->nsprings, sizeof(struct point_spring));
	struct point_spring* sprs = sprst->springs;
	
	int id = 0;
	char name[8];
	int  aid, naspr;
	double fkspr, X0, Y0, Z0;
	
	while (id < sprst->nsprings) {
		c = fscanf(fp, "%i %lf %lf %lf %lf", &naspr, &fkspr, &X0, &Y0, &Z0);
		if (c != 5) {
			ERR_MSG("Wrong spring file format\n");
		}
		
		sprs[id].naspr = naspr;
		sprs[id].fkspr = fkspr;
		sprs[id].X0 = X0;
		sprs[id].Y0 = Y0;
		sprs[id].Z0 = Z0;
		sprs[id].laspr = MYCALLOC(sprs[id].laspr, naspr, sizeof(int));
		sprs[id].ag = ag;
		
	    for (int i = 0; i < naspr; i++) {
	    	c = fscanf(fp, "%i %s", &aid, name);
	    	if (c != 2) {
				ERR_MSG("Wrong spring file format\n");
			}
			
			if (strstr(ag->atom_name[aid - 1], name) == NULL) {
				ERR_MSG("Inconsistent numbering in file %s\n"
				        "Provided atom id %i has name %s instead of %s\n", 
				         sfile, aid, ag->atom_name[aid - 1], name);
			}
			
			sprs[id].laspr[i] = aid - 1;
	    }
		
		id++;
	}
	
	fclose(fp);
	return sprst;
}

void free_springset_pairs(struct springset_pairs** sprst) 
{
	if (*sprst != NULL) {
		free((*sprst)->springs);
		free(*sprst);
		*sprst = NULL;
	}
}

void free_springset_points(struct springset_points** sprst) 
{
	if (*sprst != NULL) {
		for (int i = 0; i < (*sprst)->nsprings; i++) { 
			free((*sprst)->springs[i].laspr);
		}
		free((*sprst)->springs);
		free(*sprst);
		*sprst = NULL;
	}
}

void springeng_point(struct springset_points *sprst, double* een)
{
	int    i, i1, i2, nat;
	double xtot, ytot, ztot, fk;
	struct mol_vector3 g;
	struct mol_atom_group *ag;

	for (i = 0; i < sprst->nsprings; i++)
	{
		nat = sprst->springs[i].naspr;
		
		if(nat > 0)
		{
			xtot = 0.0;
			ytot = 0.0;
			ztot = 0.0;
			ag = sprst->springs[i].ag;

			for (i1 = 0; i1 < nat; i1++)
			{
				i2 = sprst->springs[i].laspr[i1];
				xtot += ag->coords[i2].X;
				ytot += ag->coords[i2].Y;
				ztot += ag->coords[i2].Z;
			}

			xtot = xtot / nat - sprst->springs[i].X0;
			ytot = ytot / nat - sprst->springs[i].Y0;
			ztot = ztot / nat - sprst->springs[i].Z0;

			fk = sprst->springs[i].fkspr;
			(*een) += fk * (xtot * xtot + ytot * ytot + ztot * ztot);

			fk = 2 * fk / nat;
			g.X = xtot*fk;
			g.Y = ytot*fk;
			g.Z = ztot*fk;
			
			for (i1 = 0; i1 < nat; i1++)
			{
				i2 = sprst->springs[i].laspr[i1];
				MOL_VEC_SUB(ag->gradients[i2], ag->gradients[i2], g);
			}
		}
	}
}

void springeng_pair(struct springset_pairs *sprst, double* een)
{
	int    i, i1, i2;
	double xtot, ytot, ztot, fk, d, d2, ln, er, coef, delta;
	struct mol_vector3 g;
	struct mol_atom_group *ag;

	for (i = 0; i < sprst->nsprings; i++)
	{
		ag = sprst->springs[i].ag;

		ln  = sprst->springs[i].lnspr;
		er  = sprst->springs[i].erspr;
		fk  = sprst->springs[i].fkspr / 2.0;
		
		i1 = sprst->springs[i].laspr[0];
		i2 = sprst->springs[i].laspr[1];
		
		xtot = ag->coords[i2].X - ag->coords[i1].X;
		ytot = ag->coords[i2].Y - ag->coords[i1].Y;
		ztot = ag->coords[i2].Z - ag->coords[i1].Z;
		
		d2 = xtot*xtot + ytot*ytot + ztot*ztot;
		d  = sqrt(d2);

		delta = fabs(d - ln);
		delta = (delta > er) ? ((delta - er) * delta / (d - ln)) : 0.0;
		
		//(*een) += fk * (d - ln) * (d - ln);
		(*een) += fk * delta * delta;
		//coef = fk * 2 * (1.0 - ln / d);
		coef = fk * 2.0 * delta / d;

		g.X = - coef * xtot;
		g.Y = - coef * ytot;
		g.Z = - coef * ztot;
		
		MOL_VEC_SUB(ag->gradients[i1], ag->gradients[i1], g);
		MOL_VEC_ADD(ag->gradients[i2], ag->gradients[i2], g);
	}
}

void usage_message(char** argv) {
	printf("\nUsage %s [ prm ] [ rtf ] [ -pdb ] [ -psf ]\n"
	         "         [ -rec-pdb ] [ -rec-psf  ] [ -lig-pdb  ] [ -lig-psf ]\n"
	         "         [ -lig-jsn ] [ -protocol ] [ -gbsa     ] [ -verbose ]\n"
	         "         [ -help    ] [ -nsteps   ] [ -out      ]\n\n", argv[0]);
	
	help_message();
}

void help_message(void) {
	printf("--prm      - Parameter file (required)\n"
	       "--rtf      - Topology file (required)\n"
	       "--pdb      - Full molecule pdb file\n"
	       "--psf      - Full molecule psf file\n"
	       "--rec-pdb  - Receptor pdb file\n"
	       "--rec-psf  - Receptor psf file\n"
	       "--lig-pdb  - Ligand pdb file\n"
	       "--lig-psf  - Ligand psf file\n"
	       "--lig-jsn  - Ligand json file\n"
	       "--out      - Minimized structure (default: min.pdb)\n"
	       "--protocol - Protocol file\n"
	       "--nsteps   - Number of minimization steps (default: 1000)\n"
	       "--gbsa     - Turn GBSA on (default: off)\n"
	       //"--verbose  - Verbosity\n"
	       "--help     - This message\n\n"
	       
	       "Prm and rft are required. One of the following must be provided:\n\n"
	       "\t * pdb and psf of the full molecule (--pdb --psf)\n"
	       "\t * OR rec-pdb, rec-psf, lig-pdb and lig-psf\n"
	       "\t * OR rec-pdb, rec-psf, lig-jsn\n\n"
	       
	       "Protocol file format:\n\n"
	       "\tnumber_of_steps1 fixed_atoms_file1 pair_springs_file1 point_springs_file1\n"
	       "\tnumber_of_steps2 fixed_atoms_file2 pair_springs_file2 point_springs_file2\n"
	       "\t...\n\n"
	       
	       "Pairwise springs file format:\n\n"
	       "\tnumber_of_springs\n"
	       "\tlength relative_error force_constant atom1_id(PDB) atom1_name atom2_id(PDB) atom2_name\n"
	       "\t...\n\n"
	       
	       "Point springs (attached to a single point) file format:\n\n"
	       "\tnumber_of_springs\n"
	       "\tnumber_of_atoms_attached force_constant X0 Y0 Z0\n\tatom1_id(from PDB) atom1_name\n\tatom2_id(from PDB) atom2_name\n"
	       "\t...\n\n"
	       
	       "Fixed atoms file format: simply a slice of pdb containing fixed atoms\n\n"
	       
	       "If you wish to skip any of the files in the protocol file,\n"
	       "just replace it with a dot. For instance if no fixed atoms\n"
	       "are needed then the line should be \n"
	       "\t1000 . spring_pairs spring_points\n\n"
	       
	       "CAREFUL: If ligand and receptor are specified separately,\n"
	       "\tthen ligand atom ids need to be increased by the number of receptor atoms.\n"
	       "\tCorrectness is verified by comparing atom names.\n\n");   
}

