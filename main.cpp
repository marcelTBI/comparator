#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <string>
#include <map>
#include <vector>
#include <algorithm>

extern "C" {
  #include "comparator_cmdline.h"
}

using namespace std;

#define LMINBASE 500

// reads a line no matter how long
char* my_getline(FILE *fp)
{
  char s[512], *line, *cp;
  line = NULL;
  do {
    if(fgets(s, 512, fp) == NULL) break;
    cp = strchr(s, '\n');
    if(cp != NULL) *cp = '\0';
    if(line == NULL) line = (char *) calloc(strlen(s) + 1, sizeof(char));
    else line = (char *) realloc(line, strlen(s) + strlen(line) + 1);
    strcat (line, s);
  } while (cp == NULL);
  return (line);
}

struct locstruct {
  int num;
  float energy;
  int many;         // how many walks ended here

  bool operator<(const locstruct &second) const {
    return num<second.num;
  }

} _locstruct;


// parse loc file!
int Parsefile (FILE *fp, map<string, locstruct > &mapping)
{
  char *line = NULL, *p, sep[] = " \t";
  int count = 0;

  line = my_getline(fp);
  free(line); line=NULL;

  for (count = 0, line = my_getline(fp); line != NULL; count++, line = my_getline(fp)) {
    p = strtok(line, sep);
    locstruct l;
    sscanf(p, "%d",  &l.num);

    p = strtok(NULL, sep);
    string str = p;
    p = strtok(NULL, sep);
    sscanf(p, "%f",  &l.energy);
    p = strtok(NULL, sep);
    if (p!=NULL) sscanf(p, "%d",  &l.many);
    else l.many = -1;
    mapping.insert(make_pair(str, l));
    if (line != NULL) free(line);
  }
  if (line != NULL) free(line);

  return count;
}

int main(int argc, char **argv)
{
  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "Argument parsing problem.");
    exit(EXIT_FAILURE);
  }

  FILE *barrier;
  FILE *rnaloc;
  FILE *output;

  bool barrier_file;


  // standard or barrier?
  if (args_info.standard_given) {
    barrier_file = false;
    barrier = fopen(args_info.standard_arg, "r");
    if (barrier==NULL) {
      fprintf(stderr, "Cannot open file %s\n", args_info.standard_arg);
      exit(EXIT_FAILURE);
    }
  } else {
    barrier = fopen(args_info.barriers_arg, "r");
    barrier_file = true;
    if (barrier==NULL) {
      fprintf(stderr, "Cannot open file %s\n", args_info.barriers_arg);
      exit(EXIT_FAILURE);
    }
  }

  rnaloc = fopen(args_info.rnalocmin_arg, "r");
  if (rnaloc==NULL) {
    fprintf(stderr, "Cannot open file %s\n", args_info.rnalocmin_arg);
    fclose(barrier);
    exit(EXIT_FAILURE);
  }

  if (args_info.output_given) {
    output = fopen(args_info.output_arg, "w");
    if (output==NULL) {
      fprintf(stderr, "Cannot open file %s\n", args_info.output_arg);
      fclose(barrier);
      fclose(rnaloc);
      exit(EXIT_FAILURE);
    }
  } else {
    output = stdout;
  }

  map<string, locstruct > map_loc;

  int first_nfound = -1;
  int num_found = 0;

  // parse loc file
  int loc_mins = Parsefile(rnaloc, map_loc);

  // open chain structure
  string str_open = "";
  int len = map_loc.begin()->first.length();
  for (int i=0; i<len; i++) str_open += ".";
  int ret = 0;

  char *line = NULL, *p, sep[] = " ";
  int count = 0;

  // work with barriers
    // sequence dump
  line = my_getline(barrier);
  fprintf(output, "N.(barriers) N.(locmin) %s energy    Bsize FatherBS Grad.BS    GradBBS #Walks  BS/W GBS/W\n", line+5);
  free(line); line=NULL;

  // basin statistics for missing LM
  int BSmiss = 0;
  int FBSmiss = 0;
  int GRBSmiss = 0;
  double FGRmiss = 0;
  int BScountmiss = 0;

  // also basin statistics for barr
  int BSall = 0;
  int GRBSall = 0;
  int FBSall = 0;
  double FGRall = 0;

  // also minimal energy
  float min_energy = 1e10;

  // structures
  for (count = 0, line = my_getline(barrier); line != NULL; count++, line = my_getline(barrier)) {
    p = strtok(line, sep);
    int number;
    sscanf(p, "%d",  &number);

    p = strtok(NULL, sep);
    string str = p;
    p = strtok(NULL, sep);
    float energy;
    sscanf(p, "%f", &energy);

    // apply energy range
    if (count > 0 && energy > min_energy + args_info.erange_arg) {
      break;
    }

    min_energy = min(min_energy, energy);

    int father;
    float en_diff;
    int bsize=0;          /* # of structures in this lmin  == basin size */
    //double eff_bsize;      /* bsize without children */
    int fathers_bsize=0;  /* bsize of father at merging-time */
    double F;              /* free energy of lmin i and its children */
    //double F_eff;          /* free engery of lmin (alone) */
    //double Z;              /* partition function of lmin */
    //double Z_eff;          /* Z of lmin (alone) */
    int Gr_bsize=0;       /* gradient basin size */
    double FGr;            /* F of gradient basin */

    if (barrier_file) {
      p = strtok(NULL, sep); sscanf(p, "%d", &father);
      p = strtok(NULL, sep); sscanf(p, "%f", &en_diff);
      p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%d", &bsize);
      p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%d", &fathers_bsize);
      p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%lf", &F);
      p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%d", &Gr_bsize);
      p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%lf", &FGr);}}}}}
    } else {
      p = strtok(NULL, sep); sscanf(p, "%d", &Gr_bsize);
    }

    double _kT = 0.00198717*(273.15 + 37);
    //printf("%5.4f", FGr);
    FGr = exp((min_energy-FGr)/_kT)*Gr_bsize;

    map<string, locstruct >::iterator it;
    it=map_loc.find(str);

    // open_chain?
    if ((args_info.open_barr_flag || args_info.open_loc_flag) && str==str_open) {
      if (args_info.open_barr_flag) {
        ret = number;
      } else {
        if (it!=map_loc.end()) {
          ret = it->second.num;
        }
      }
    }

    // print output
    if (it != map_loc.end()) {
      num_found ++;
      if (it->second.energy != energy) {
        fprintf(stderr, "WARNING: energies does not agree (%f != %f @ %d)\n", it->second.energy, energy, count);
      }
      loc_mins--;
      if (Gr_bsize==0) {
        fprintf(output, "%12d %10d %s %6.2f\n",
                number, it->second.num, str.c_str(), energy);
      } else {
        fprintf(output, "%12d %10d %s %6.2f %8d %8d %7d %10.1f %6d %5.1f %5.1f\n",
                number, it->second.num, str.c_str(), energy, bsize, fathers_bsize, Gr_bsize, FGr, it->second.many, bsize/(float)it->second.many, Gr_bsize/(float)it->second.many);
      }
      map_loc.erase(it);
    } else {
      if (first_nfound<0) first_nfound = count;
      if (Gr_bsize==0) {
        fprintf(output, "%12d            %s %6.2f\n", number, str.c_str(), energy);
      } else {
        BSmiss += bsize;
        FBSmiss += fathers_bsize;
        GRBSmiss += Gr_bsize;
        FGRmiss += FGr;
        BScountmiss++;
        fprintf(output, "%12d            %s %6.2f %8d %8d %7d %10.1f\n", number, str.c_str(), energy, bsize, fathers_bsize, Gr_bsize, FGr);
      }
    }

    // statistics about all
    if (Gr_bsize!=0) {
      BSall+=bsize;
      GRBSall+=Gr_bsize;
      FBSall+=fathers_bsize;
      FGRall+=FGr;
    }

    if (line != NULL) free(line);
  }
  if (line != NULL) free(line);

  // display basin statistics
  if (GRBSall!= 0) {
    for (int i=0; i<13+13+len+8-24-4;i++) fprintf(output, " ");
    fprintf(output, "%3d total missing basins: %8d %8d %7d %10.1f\n", BScountmiss, BSmiss, FBSmiss, GRBSmiss, FGRmiss);
    for (int i=0; i<13+13+len+8-24-4;i++) fprintf(output, " ");
    fprintf(output, "  average missing basins: %8.2lf %8.2lf %7.2lf %10.1f\n", BSmiss/(double)BScountmiss, FBSmiss/(double)BScountmiss, GRBSmiss/(double)BScountmiss, FGRmiss/(double)BScountmiss);
    for (int i=0; i<13+13+len+8-24-4;i++) fprintf(output, " ");
    fprintf(output, "  total           basins: %8d %8d %7d %10.1f\n", BSall, FBSall, GRBSall, FGRall);
    for (int i=0; i<13+13+len+8-24-4;i++) fprintf(output, " ");
    fprintf(output, "  percent missing basins: %8.2lf%% %7.2lf%% %6.2lf%%    %6.2f%%\n", BSmiss/(double)BSall*100.0, FBSmiss/(double)FBSall*100.0, GRBSmiss/(double)GRBSall*100.0, FGRmiss/FGRall*100.0);
  }

  if (loc_mins!=0) {
    fprintf(output, "\n");
    //fprintf(stderr, "WARNING: %d local minima have not been found by barriers\n", loc_mins);

    // sort 'em in vector and print
    vector<pair<locstruct, string> > out;
    for (map<string, locstruct >::iterator it=map_loc.begin(); it!=map_loc.end(); it++) {
      out.push_back(make_pair(it->second, it->first));
      //fprintf(output, "             %10d %s %6.2f\n", it->second.first, it->first.c_str(), it->second.second);
    }
    sort(out.begin(), out.end());
    for (unsigned int i=0; i<out.size(); i++) {
      if (out[i].first.energy > min_energy + args_info.erange_arg) break;
      fprintf(output, "             %10d %s %6.2f %8d\n", out[i].first.num, out[i].second.c_str(), out[i].first.energy, out[i].first.many);
    }
  }

  // open loc?
  if (args_info.open_loc_flag && ret==0) {
    if (map_loc.find(str_open)!=map_loc.end()) {
      ret = map_loc[str_open].num;
    } else {
      fprintf(stderr, "WARNING: RNAlocmin haven\'t found open chain!\n");
    }
  }

  // open barr
  if (args_info.open_barr_flag && ret==0) {
      fprintf(stderr, "WARNING: barriers haven\'t found open chain!\n");
  }

  // print 6 numbers: number of first not found minimum + num of found minima + percent of missed Gr basins, mfe, missed fgr, number of LM
  printf("%3d %3d %5.2f %7.2f %5.2f %3d\n", first_nfound+1, num_found, GRBSmiss/(double)GRBSall*100.0, min_energy, FGRmiss/FGRall*100.0, count);

  // output percent of missing basins to another file
  if (args_info.output_basin_given) {
    FILE *fil = fopen(args_info.output_basin_arg, "w");

    fprintf(fil, "%.2f\n", GRBSmiss/(double)GRBSall*100.0);

    fclose(fil);
  }

  fclose(output);
  fclose(barrier);
  fclose(rnaloc);
  cmdline_parser_free(&args_info);

  if (args_info.open_loc_flag || args_info.open_barr_flag) return ret;
  else return loc_mins;
}
