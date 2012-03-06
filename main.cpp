#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

} _locstruct;

int Parsefile (FILE *fp, map<string, locstruct > &mapping)
{
  char *line = NULL, *p, sep[] = " ";
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
    sscanf(p, "%d",  &l.many);

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

  barrier = fopen(args_info.barriers_arg, "r");
  if (barrier==NULL) {
    fprintf(stderr, "Cannot open file %s\n", args_info.barriers_arg);
    exit(EXIT_FAILURE);
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
  int last_found;

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
  fprintf(output, "N.(barriers) N.(locmin) %s energy    Bsize FatherBS Grad.BS #Walks  BS/W GBS/W\n", line+5);
  free(line); line=NULL;

  // basin statistics
  int BSall = 0;
  int FBSall = 0;
  int GRBSall = 0;
  int BScount = 0;

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

    p = strtok(NULL, sep); sscanf(p, "%d", &father);
    p = strtok(NULL, sep); sscanf(p, "%f", &en_diff);
    p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%d", &bsize);
    p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%d", &fathers_bsize);
    p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%lf", &F);
    p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%d", &Gr_bsize);
    p = strtok(NULL, sep); if (p!=NULL) {sscanf(p, "%lf", &FGr);}}}}}

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
      last_found = it->second.num;
      if (it->second.energy != energy) {
        fprintf(stderr, "WARNING: energies does not agree (%f != %f @ %d)\n", it->second.energy, energy, count);
      }
      loc_mins--;
      if (bsize==0) {
        fprintf(output, "%12d %10d %s %6.2f\n",
                number, it->second.num, str.c_str(), energy);
      } else {
        fprintf(output, "%12d %10d %s %6.2f %8d %8d %7d %6d %5.1f %5.1f\n",
                number, it->second.num, str.c_str(), energy, bsize, fathers_bsize, Gr_bsize, it->second.many, bsize/(float)it->second.many, Gr_bsize/(float)it->second.many);
      }
      map_loc.erase(it);
    } else {
      if (first_nfound<0) first_nfound = count;
      if (bsize==0) {
        fprintf(output, "%12d            %s %6.2f\n", number, str.c_str(), energy);
      } else {
        BSall += bsize;
        FBSall += fathers_bsize;
        GRBSall += Gr_bsize;
        BScount++;
        fprintf(output, "%12d            %s %6.2f %8d %8d %7d\n", number, str.c_str(), energy, bsize, fathers_bsize, Gr_bsize);
      }
    }

    if (line != NULL) free(line);
  }
  if (line != NULL) free(line);

  // display basin statistics
  if (BSall != 0) {
    for (int i=0; i<13+13+len+8-24-4;i++) fprintf(output, " ");
    fprintf(output, "%3d total missing basins: %8d %8d %7d\n", BScount, BSall, FBSall, GRBSall);
    for (int i=0; i<13+13+len+8-24-4;i++) fprintf(output, " ");
    fprintf(output, "  average missing basins: %8.2lf %8.2lf %7.2lf\n", BSall/(double)BScount, FBSall/(double)BScount, GRBSall/(double)BScount);
  }

  if (loc_mins!=0) {
    fprintf(output, "\n");
    //fprintf(stderr, "WARNING: %d local minima have not been found by barriers\n", loc_mins);

    // sort 'em in vector and print
    vector<pair<int, string> > out;
    for (map<string, locstruct >::iterator it=map_loc.begin(); it!=map_loc.end(); it++) {
      out.push_back(make_pair(it->second.num, it->first));
      //fprintf(output, "             %10d %s %6.2f\n", it->second.first, it->first.c_str(), it->second.second);
    }
    sort(out.begin(), out.end());
    for (unsigned int i=0; i<out.size(); i++) {
      fprintf(output, "             %10d %s %6.2f\n", out[i].first, out[i].second.c_str(), map_loc[out[i].second].energy);
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

  // print two numbers: num of found minima + number of first not found minimum
  printf("%3d %3d\n", first_nfound+1, last_found);

  fclose(output);
  fclose(barrier);
  fclose(rnaloc);
  cmdline_parser_free(&args_info);

  if (args_info.open_loc_flag || args_info.open_barr_flag) return ret;
  else return loc_mins;
}
