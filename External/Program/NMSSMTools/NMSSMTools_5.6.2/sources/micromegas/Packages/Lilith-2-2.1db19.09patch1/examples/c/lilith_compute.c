#include <Python.h>
#include "lilith.h"

int main(int argc, char* argv[])
{
    // Initializing Python/C interface
    Py_Initialize();

    // Parameters:
    // * Experimental list path
    // * Path to a user input file
    char experimental_input[] = "../../data/latest.list";
//    char XMLinputpath[] = "userinput/example_couplings.xml";
  
    // Accessible outputs for a given point
    // * Output of the reduced couplings
    // * Output of the various contributions to the total -2LogL in a XML format
    // * Output of the total -2LogL in a SLHA format
    // * Output of the signal strengths
    char output_couplings[] = "lilith_couplings_output.xml";
    char output_XML[] = "lilith_likelihood_output.xml";
    char output_SLHA[] = "lilith_likelihood_output.slha";
    char output_mu[] = "lilith_mu_output.xml";
  
  
  
    // Creating an object of the class Lilith: lilithcalc
    PyObject* lilithcalc = initialize_lilith(experimental_input);
  
  
    // Creating an XML input string for the reduced coupling mode
    // signalstrength mode would work as well
    char XMLinputstring[6000]="";
    char buffer[100];
  
    double mh = 125.;
    double CU = 1.;
    double CD = 1.;
    double CV = 1.;
    double CGa = 1.;
    double Cg = 1.;
    double BRinv = 0.;
    double BRund = 0.;
    char precision[] = "BEST-QCD";

    sprintf(buffer,"<?xml version=\"1.0\"?>\n");
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<lilithinput>\n");
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<reducedcouplings>\n");
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<mass>%f</mass>\n", mh);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"tt\">%f</C>\n", CU);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"cc\">%f</C>\n", CU);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"bb\">%f</C>\n", CD);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"tautau\">%f</C>\n", CD);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"ZZ\">%f</C>\n", CV);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"WW\">%f</C>\n", CV);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"gammagamma\">%f</C>\n", CGa);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<C to=\"gg\">%f</C>\n", Cg);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<extraBR>\n");
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<BR to=\"invisible\">%f</BR>\n", BRinv);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<BR to=\"undetected\">%f</BR>\n", BRund);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"</extraBR>\n");
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"<precision>%s</precision>\n",precision);
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"</reducedcouplings>\n");
    strcat(XMLinputstring, buffer);
    sprintf(buffer,"</lilithinput>\n");
    strcat(XMLinputstring, buffer);

    // Reading user input XML string
    lilith_readuserinput(lilithcalc, XMLinputstring);

    // Getting -2LogL
    float my_likelihood;
    my_likelihood = lilith_computelikelihood(lilithcalc);
    printf("-2*log(L) = %lf\n", my_likelihood);
  
  
    // Getting exp_ndf
    int exp_ndf;
    exp_ndf = lilith_exp_ndf(lilithcalc);
    printf("exp_ndf = %i\n", exp_ndf);
  
    // Writing the likelihood results in XML file
    lilith_likelihood_output(lilithcalc, output_XML, 0);
  
    // Writing the likelihood results in SLHA file
    lilith_likelihood_output(lilithcalc, output_SLHA, 1);
  
    // Writing the signal strength results in a XML file
    lilith_mu_output(lilithcalc, output_mu, 0);
  
    // Writing the couplings in a XML file
    lilith_couplings_output(lilithcalc, output_couplings);

    // Exiting Python/C interface
    Py_Finalize();
    return 0;
}

