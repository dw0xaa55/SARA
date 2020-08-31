/*
 * SARA - Spectral Analyzer for Radio Astronomy
 *
 * Description :
 * This program is meant for use with data obtained by the SALSA Telescope located in Onsala, Sweden.
 * The spectral data can be cleaned from background noise and normalized for displaying coordinates
 * of HI clouds which are calculated from the available data.
 * This program is currently in alpha state so expect it to behave buggy from time to time.
 *
 * Disclaimer  :
 * I am by no means a professional coder in C or the GTK Framework. Basically that means, I tought it
 * to myself. So expect this code not to be very clean or elegant. I plan to optimize the source code
 * in the future, which will take some time, because therefore many things have to be learned to be
 * achieved.
 * In the mean time: Have fun using this version of SARA :>
 *
 * Version     : 0.6 [alpha]
 * Author      : C. Huffenbach
 * E-Mail      : nphard1234@gmail.com 
 * Date        : 2020-06-01
 *
 */

// TODO:
/*
 * - [X] make colors for map rendering dependent on the intensity (look at mapData[][3], thats the intensity (why isnt it displayed in draw function)(maxValue of intensity = 140))
 * - [ ] maybe generate rainbow colors from new colorOfMap variable
 * - [ ] map view add zoom buttons which control zoomModifier;
 * - [ ] map view add color buttons which contoll colorModifier;
 * - [ ] clear map generator data before loading new data
 * - [ ] delete paths from monitoredFilePaths by clicking "-" sign and clear menu entry (maybe not neccessary because single fitted file for map generator) -> clean up unused monitoredFilePaths
 * - [ ] find another way to validate rMinus and rPlus/ possibly find another visualization algorithm for reading all four quadrants (if possible)
 * - [ ] add status bar with quick tips
 * - [ ] up/down button for file selector
 * - [ ] delete file->save
 * - [ ] generate and save an image of generated map
 * - [ ] about und version dialogs
 * - [ ] clean up code
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>
#include <gtk/gtk.h>
#include <gtk/gtkx.h>
#include <math.h>
#include <ctype.h>
#include <sys/mman.h>

// declare GtkWidgets here
GtkBuilder        *builder;
GtkWidget         *window;
GtkWidget         *mapWindow;
GtkWidget         *curveEditorWindow;
GtkWidget         *gaussWindow;
GtkWidget         *fixedGridGauss;
GtkWidget         *fixedGridEditor;
GtkWidget         *fixedGrid;
GtkWidget         *menubar;
GtkWidget         *clear; 
GtkWidget         *openFile;
GtkWidget         *openStack;
GtkWidget         *save;
GtkWidget         *exitProg;
GtkWidget         *blackWhite;
GtkWidget         *temperature;
GtkWidget         *custom;
GtkWidget         *curve;
GtkWidget         *map;
GtkWidget         *version;
GtkWidget         *help;
GtkWidget         *scrolledGaussValues;
GtkWidget         *scrolledWindow;
GtkWidget         *viewportGauss;
GtkWidget         *viewportTextview;
GtkWidget         *textview;
GtkWidget         *btn_deleteGaussCoords;
GtkWidget         *btn_deleteAllGaussCoords;
GtkWidget         *btn_exportGaussCoords;
GtkWidget         *btn_fileAdd;
GtkWidget         *btn_fileDel;
GtkWidget         *btn_filePreview;
GtkWidget         *btn_baseline;
GtkWidget         *btn_substractBaseline;
GtkWidget         *btn_gaussFit;
GtkWidget         *btn_saveGauss;
GtkWidget         *btn_clearCurve;
GtkWidget         *graphDraw;
GtkWidget         *btn_scaleIncrease;
GtkWidget         *btn_scaleDecrease;
GtkWidget         *lbl_Temp1;
GtkWidget         *lbl_Temp2;
GtkWidget         *lbl_Temp3;
GtkWidget         *lbl_Temp4;
GtkWidget         *lbl_Temp5;
GtkWidget         *lbl_Temp6;
GtkWidget         *lbl_Temp7;
GtkWidget         *lbl_Temp8;
GtkWidget         *lbl_Temp9;
GtkWidget         *mapDrawingArea;
GtkWidget         *scrolledDistanceValues;
GtkWidget         *viewportMap;
GtkWidget         *scrolledFileList;
GtkWidget         *viewportFileList;
GtkWidget         *lbl_header;
GtkCellRenderer   *crGaussLon;
GtkCellRenderer   *crGaussLat;
GtkCellRenderer   *crGaussVel;
GtkCellRenderer   *crGaussTemp;
GtkCellRenderer   *crGalLon;
GtkCellRenderer   *crVelObs;
GtkCellRenderer   *crDistR;
GtkCellRenderer   *crRplus;
GtkCellRenderer   *crRminus;
GtkCellRenderer   *crDistSol;
GtkCellRenderer   *crX;
GtkCellRenderer   *crY;
GtkCellRenderer   *filePath;
GtkTreeViewColumn *colGaussLon;
GtkTreeViewColumn *colGaussLat;
GtkTreeViewColumn *colGaussVel;
GtkTreeViewColumn *colGaussTemp;
GtkTreeViewColumn *filenames;
GtkTreeViewColumn *colGalLon;
GtkTreeViewColumn *colVelObs;
GtkTreeViewColumn *colDistR;
GtkTreeViewColumn *colRplus;
GtkTreeViewColumn *colRminus;
GtkTreeViewColumn *colDistSol;
GtkTreeViewColumn *colX;
GtkTreeViewColumn *colY;
GtkTextBuffer     *textBuffer;
GtkTreeIter       iter;
GtkTreeIter       mapIter;
GtkTreeIter       gaussIter;
GtkTreeModel      *modelGauss;
GtkTreeModel      *model;
GtkTreeStore      *treeStoreGauss;
GtkTreeStore      *treeStore;
GtkTreeStore      *storeMap;
GtkTreeView       *treeViewGauss;
GtkTreeView       *treeViewFiles;
GtkTreeView       *treeMapValues;
GtkTreeSelection  *gaussSelection;
GtkTreeSelection  *fileSelection;
GtkTreeSelection  *mapSelection;

// GtkWidget *btn_prevSelect;
// GtkWidget *btn_nextSelect;

// global variables
char    filename[512],multifilename[512],strLblHeader[512],monitoredFilePaths[512][512];
float   graphX[512],graphY[512],galLonData[512],velObs[512][256],temperatureObs[512][256],cloudPositionX[92160],cloudPositionY[92160],gLon;
guint   width, height;
int     scaleModifier = 2,fileListLength = 0, clickCounter=0;
Bool    previewActive = FALSE,reduced = FALSE,gaussFit=FALSE;
gdouble baseX[4], baseY[4];
double  XMatrix[4][4], YMatrix[4], XMatrixInv[4][4], baselineMatrix[4],baselineCoords[255],baselineCoordsCopy[255],gaussBellCurve[255];

int main(int argc, char *argv[]){
  // initiate gtk and default stuff
  gtk_init(&argc,&argv);

  builder           = gtk_builder_new_from_file("/usr/share/sara/sara.glade");
  window            = GTK_WIDGET (gtk_builder_get_object (builder, "window"            ));
  mapWindow         = GTK_WIDGET (gtk_builder_get_object (builder, "mapWindow"         ));
  curveEditorWindow = GTK_WIDGET (gtk_builder_get_object (builder, "curveEditorWindow" ));
  gaussWindow       = GTK_WIDGET (gtk_builder_get_object (builder, "gaussWindow"       ));
  
  g_signal_connect(window,            "destroy",            G_CALLBACK(gtk_main_quit),        NULL);

  gtk_builder_connect_signals(builder,NULL);
  
  // initiate GtkWidgets here
  fixedGrid                = GTK_WIDGET           (gtk_builder_get_object(builder, "fixedGrid"                ));
  fixedGridEditor          = GTK_WIDGET           (gtk_builder_get_object(builder, "fixedGridEditor"          ));
  fixedGridGauss           = GTK_WIDGET           (gtk_builder_get_object(builder, "fixedGridGauss"           ));
  menubar                  = GTK_WIDGET           (gtk_builder_get_object(builder, "menubar"                  ));
  clear                    = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_clear"                 ));
  openFile                 = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_open"                  ));
  openStack                = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_openStack"             ));
  save                     = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_save"                  ));
  exitProg                 = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_exit"                  ));
  blackWhite               = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_blackWhite"            ));
  temperature              = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_temperature"           ));
  custom                   = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_custom"                ));
  curve                    = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_curve"                 ));
  map                      = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_map"                   ));
  version                  = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_version"               ));
  help                     = GTK_WIDGET           (gtk_builder_get_object(builder, "mi_help"                  ));
  scrolledWindow           = GTK_WIDGET           (gtk_builder_get_object(builder, "scrolledWindow"           ));
  viewportTextview         = GTK_WIDGET           (gtk_builder_get_object(builder, "viewportTextview"         ));
  viewportGauss            = GTK_WIDGET           (gtk_builder_get_object(builder, "viewportGauss"            ));
  textview                 = GTK_WIDGET           (gtk_builder_get_object(builder, "displayTextview"          ));
  lbl_header               = GTK_WIDGET           (gtk_builder_get_object(builder, "lbl_header"               ));
  scrolledFileList         = GTK_WIDGET           (gtk_builder_get_object(builder, "scrolledFileList"         ));
  viewportFileList         = GTK_WIDGET           (gtk_builder_get_object(builder, "viewportFileList"         ));
  treeStore                = GTK_TREE_STORE       (gtk_builder_get_object(builder, "treeStore"                ));
  treeViewFiles            = GTK_TREE_VIEW        (gtk_builder_get_object(builder, "treeViewFiles"            ));
  treeViewGauss            = GTK_TREE_VIEW        (gtk_builder_get_object(builder, "treeViewGauss"            ));
  filenames                = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "filenames"                ));
  filePath                 = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "filePath"                 ));
  fileSelection            = GTK_TREE_SELECTION   (gtk_builder_get_object(builder, "fileSelection"            ));
  gaussSelection           = GTK_TREE_SELECTION   (gtk_builder_get_object(builder, "gaussSelection"           ));
  btn_fileAdd              = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_fileAdd"              ));
  btn_fileDel              = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_fileDel"              ));
  btn_filePreview          = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_filePreview"          ));
  btn_baseline             = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_baseline"             ));
  btn_substractBaseline    = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_substractBaseline"    ));
  btn_gaussFit             = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_gaussFit"             ));
  btn_saveGauss            = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_saveGauss"            ));
  btn_clearCurve           = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_clearCurve"           ));
  btn_deleteGaussCoords    = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_deleteGaussCoords"    ));
  btn_deleteAllGaussCoords = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_deleteAllGaussCoords" ));
  btn_exportGaussCoords    = GTK_WIDGET           (gtk_builder_get_object(builder, "btn_exportGaussCoords"    ));
  graphDraw                = GTK_WIDGET           (gtk_builder_get_object(builder, "graphDraw"                ));
  btn_scaleIncrease        = GTK_WIDGET           (gtk_builder_get_object(builder, "btnTempIncrease"          ));
  btn_scaleDecrease        = GTK_WIDGET           (gtk_builder_get_object(builder, "btnTempDecrease"          ));
  lbl_Temp1                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp1"                 ));
  lbl_Temp2                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp2"                 ));
  lbl_Temp3                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp3"                 ));
  lbl_Temp4                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp4"                 ));
  lbl_Temp5                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp5"                 ));
  lbl_Temp6                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp6"                 ));
  lbl_Temp7                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp7"                 ));
  lbl_Temp8                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp8"                 ));
  lbl_Temp9                = GTK_WIDGET           (gtk_builder_get_object(builder, "lblTemp9"                 ));
  mapDrawingArea           = GTK_WIDGET           (gtk_builder_get_object(builder, "mapDrawingArea"           ));
  scrolledDistanceValues   = GTK_WIDGET           (gtk_builder_get_object(builder, "scolledDistanceValues"    ));
  scrolledGaussValues      = GTK_WIDGET           (gtk_builder_get_object(builder, "scrolledGaussValues"      ));
  viewportMap              = GTK_WIDGET           (gtk_builder_get_object(builder, "viewportMap"              ));
  storeMap                 = GTK_TREE_STORE       (gtk_builder_get_object(builder, "storeMap"                 ));
  treeStoreGauss           = GTK_TREE_STORE       (gtk_builder_get_object(builder, "treeStoreGauss"           ));
  treeMapValues            = GTK_TREE_VIEW        (gtk_builder_get_object(builder, "treeMapValues"            ));
  colGalLon                = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colGalLon"                ));
  colVelObs                = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colVelObs"                ));
  colDistR                 = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colDistR"                 ));
  colRplus                 = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colRplus"                 ));
  colRminus                = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colRminus"                ));
  colDistSol               = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colDistSol"               ));
  colX                     = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colX"                     ));
  colY                     = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colY"                     ));
  colGaussLon              = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colGaussLon"              ));
  colGaussLat              = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colGaussLat"              ));
  colGaussVel              = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colGaussVel"              ));
  colGaussTemp             = GTK_TREE_VIEW_COLUMN (gtk_builder_get_object(builder, "colGaussTemp"             ));
  crGalLon                 = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crGalLon"                 ));
  crVelObs                 = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crVelObs"                 ));
  crDistR                  = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crDistR"                  ));
  crRplus                  = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crRplus"                  ));
  crRminus                 = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crRminus"                 ));
  crDistSol                = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crDistSol"                ));
  crX                      = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crX"                      ));
  crY                      = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crY"                      ));
  crGaussLon               = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crGaussLon"               ));
  crGaussLat               = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crGaussLat"               ));
  crGaussVel               = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crGaussVel"               ));
  crGaussTemp              = GTK_CELL_RENDERER    (gtk_builder_get_object(builder, "crGaussTemp"              ));

  // btn_prevSelect = GTK_WIDGET(gtk_builder_get_object(builder,"btnPrevSelection"));
  // btn_nextSelect = GTK_WIDGET(gtk_builder_get_object(builder,"btnNextSelection"));

  fileSelection          = gtk_tree_view_get_selection (GTK_TREE_VIEW( treeViewFiles  ));
  textBuffer             = gtk_text_view_get_buffer    (GTK_TEXT_VIEW( textview       ));
  mapSelection           = gtk_tree_view_get_selection (GTK_TREE_VIEW( treeMapValues  ));
  gaussSelection         = gtk_tree_view_get_selection (GTK_TREE_VIEW( treeViewGauss  ));
  
  gtk_tree_view_column_add_attribute (filenames,    filePath,    "text", 0);
  gtk_tree_view_column_add_attribute (colGalLon,    crGalLon,    "text", 0);
  gtk_tree_view_column_add_attribute (colVelObs,    crVelObs,    "text", 1);
  gtk_tree_view_column_add_attribute (colDistR,     crDistR,     "text", 2);
  gtk_tree_view_column_add_attribute (colRplus,     crRplus,     "text", 3);
  gtk_tree_view_column_add_attribute (colRminus,    crRminus,    "text", 4);
  gtk_tree_view_column_add_attribute (colDistSol,   crDistSol,   "text", 5);
  gtk_tree_view_column_add_attribute (colX,         crX,         "text", 6);
  gtk_tree_view_column_add_attribute (colY,         crY,         "text", 7);
  gtk_tree_view_column_add_attribute (colGaussLon,  crGaussLon,  "text", 0);
  gtk_tree_view_column_add_attribute (colGaussLat,  crGaussLat,  "text", 1);
  gtk_tree_view_column_add_attribute (colGaussVel,  crGaussVel,  "text", 2);
  gtk_tree_view_column_add_attribute (colGaussTemp, crGaussTemp, "text", 3);
  
  g_object_unref(builder);

  gtk_widget_set_events(graphDraw,GDK_BUTTON_PRESS_MASK);
  
  GdkColor color;
  color.red = 0xDDDD;
  color.green = 0xDDDD;
  color.blue = 0xDDDD;
  gtk_widget_modify_bg(GTK_WIDGET(window), GTK_STATE_NORMAL, &color);

  gtk_window_set_keep_above(GTK_WINDOW(window),TRUE);
  gtk_widget_show(window);
  gtk_main();
  return EXIT_SUCCESS;
}

// all syntax and logic functions go here

// Menubar Item -> Clear
void on_mi_clear_activate(GtkMenuItem *m){
  // clear all variables and textview
  reduced = FALSE;
  previewActive = FALSE;
  gtk_widget_set_sensitive(btn_scaleIncrease,FALSE);
  gtk_widget_set_sensitive(btn_scaleDecrease,FALSE);

  // gtk_widget_set_sensitive(btn_prevSelect,TRUE);
  // gtk_widget_set_sensitive(btn_nextSelect,TRUE);

  for(int i=0;i<512;i++){
    graphX[i]=0;
    graphY[i]=0;
  }
  gtk_widget_queue_draw(graphDraw);
  memset(filename,0,strlen(filename));
  memset(strLblHeader,0,strlen(strLblHeader));
  gtk_text_buffer_set_text(textBuffer,(const gchar *) "",(gint) -1);
  gtk_label_set_text((GtkLabel *)lbl_header,(const gchar *)"        --- NO DATA LOADED ---");
  gtk_tree_store_clear(treeStore);
  gtk_tree_store_clear(treeStoreGauss);
  gtk_tree_store_clear(storeMap);
  fileListLength=0;
}

void on_mi_clear_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_clear_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> open
// read file, interprete it and display it in textview

float gLat;

void loadTextFile(char *file[]){
  fileListLength=0;
  Bool headerActive = FALSE;
  int lineNumber = 0;
  gtk_text_buffer_set_text(textBuffer,(const gchar *) "",(gint) -1);
  FILE * fPointer;
  fPointer = fopen((const char *)file,"r");
  char singleLine[256];
  char temp[32];
  int half=0;

  while(!feof(fPointer)){
    Bool xToken = TRUE;
    fgets(singleLine, 256, fPointer);
    // read line and separate header from data
    // glon und glat lesen und speichern 
    if(singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'G' && singleLine[3] == 'L' && singleLine[4] == 'O' && singleLine[5] == 'N' && singleLine[6] == '='){
      for(int i=7;i<strlen(singleLine);i++){
        if(singleLine[i]=='.'){
          temp[i-7]=',';
        }else{
          temp[i-7]=singleLine[i];
        }
        gLon = strtof((const char *restrict) temp,NULL);
      }
      memset(temp,0,strlen(temp));
    }

    if(singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'G' && singleLine[3] == 'L' && singleLine[4] == 'A' && singleLine[5] == 'T' && singleLine[6] == '='){
      for(int i=7;i<strlen(singleLine);i++){
        if(singleLine[i]=='.'){
          temp[i-7]=',';
        }else{
          temp[i-7]=singleLine[i];
        }
        gLat = strtof((const char *restrict) temp,NULL);
      }
      memset(temp,0,strlen(temp));
    }

    if(singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'B' && singleLine[3] == 'E' && singleLine[4] == 'G' && singleLine[5] == 'I' && singleLine[6] == 'N'){
      headerActive = True;
    }else if(singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'D' && singleLine[3] == 'A' && singleLine[4] == 'T' && singleLine[5] == 'E' && singleLine[6] == '=' ||
             singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'G' && singleLine[3] == 'L' && singleLine[4] == 'O' && singleLine[5] == 'N' && singleLine[6] == '=' ||
             singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'G' && singleLine[3] == 'L' && singleLine[4] == 'A' && singleLine[5] == 'T' && singleLine[6] == '=' ||
             singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'D' && singleLine[3] == 'A' && singleLine[4] == 'T' && singleLine[5] == 'A' && singleLine[6] == ' '){
      strncat(strLblHeader,(const char * restrict) &singleLine,strlen((const gchar *)singleLine));
      gtk_label_set_text((GtkLabel *) lbl_header,(const gchar *) strLblHeader);
    }else if(singleLine[0] == '#' && singleLine[1] == ' ' && singleLine[2] == 'E' && singleLine[3] == 'N' && singleLine[4] == 'D'){
      headerActive = False;
    }else if(!headerActive){
      // split columns into x and y values
      for(int i=0;i<strlen(singleLine);i++){
        if(singleLine[i]=='.'){
          temp[i]=',';
        }else{
          temp[i]=singleLine[i];
        }
        char *fend;
        graphX[lineNumber]=strtof((const char *restrict)temp,&fend);
        graphY[lineNumber]=strtof((const char *restrict)fend,NULL);
      }
      memset(temp,0,strlen(temp));
      lineNumber++;
      gtk_text_buffer_insert_at_cursor(textBuffer,(const gchar *) singleLine, (gint) -1);
    }else{
      // do nothing here because the rest are comments
    }
  }
  gtk_widget_queue_draw(graphDraw);
  fclose(fPointer);
}

void on_mi_open_activate(GtkMenuItem *m){
  // open file dialog
  GtkWidget *dialog;
  dialog=gtk_file_chooser_dialog_new("Choose a spectrum file",GTK_WINDOW(window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_OK,GTK_RESPONSE_OK,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,NULL);
  gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(dialog),FALSE);
  gtk_widget_show_all(dialog);
  gint resp=gtk_dialog_run(GTK_DIALOG(dialog));
  if(resp == GTK_RESPONSE_OK){
    strcpy(filename,gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog)));
    gtk_tree_store_append(treeStore,&iter,NULL);
    gtk_tree_store_set(treeStore,&iter,0,(char **)filename,-1);
    // put filename in monitoredFilePaths
    strcpy(monitoredFilePaths[fileListLength],filename);
    fileListLength++;
  }else{
  }
  gtk_widget_destroy(dialog);
  // gtk_widget_set_sensitive(btn_prevSelect,TRUE);
  // gtk_widget_set_sensitive(btn_nextSelect,TRUE);
}

void on_mi_open_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_open_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// read file paths from GSList read from file chooser below
void displayPath(gpointer *data){
  strcpy(filename,(const char *)data);
  gtk_tree_store_append(treeStore,&iter,NULL);
  gtk_tree_store_set(treeStore,&iter,0,(char **)filename,-1);
  // put filename in monitoredFilePaths
  strcpy(monitoredFilePaths[fileListLength],filename);
  fileListLength++;
}

// Menubar Item -> open Stack
void on_mi_openStack_activate(GtkMenuItem *m){
  // open folder dialog
  GtkWidget *dialog;
  dialog=gtk_file_chooser_dialog_new("Choose spectrum files",GTK_WINDOW(window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_OK,GTK_RESPONSE_OK,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,NULL);
  gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(dialog),TRUE);
  gtk_widget_show_all(dialog);
  gint resp=gtk_dialog_run(GTK_DIALOG(dialog));
  if(resp == GTK_RESPONSE_OK){
    GSList *tempFilelist = gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(dialog));
    guint glistLen = g_slist_length(tempFilelist);
    g_slist_foreach(tempFilelist,(GFunc)displayPath,NULL);
  }else{
  }
  gtk_widget_destroy(dialog);
}

void on_mi_openStack_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_openStack_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> save
void on_mi_save_activate(GtkMenuItem *m){
  // menubar item entry was clicked
}

void on_mi_save_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_save_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> exit
void on_mi_exit_activate(GtkMenuItem *m){
  // quit program
  gtk_main_quit();
}

void on_mi_exit_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_exit_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> color black white
void on_mi_blackWhite_activate(GtkMenuItem *m){
  // menubar item entry was clicked
}

void on_mi_blackWhite_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_blackWhite_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> color temperature
void on_mi_temperature_activate(GtkMenuItem *m){
  // menubar item entry was clicked
}

void on_mi_temperature_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_temperature_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> color custom
void on_mi_custom_activate(GtkMenuItem *m){
  // menubar item entry was clicked
}

void on_mi_custom_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_custom_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> window curve
void on_mi_curve_activate(GtkMenuItem *m){
  gtk_widget_show(curveEditorWindow);
  gtk_widget_show(gaussWindow);
}

void on_mi_curve_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_curve_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> window map
double mapData[2000][4];
float colorOfMap[92160];

void on_mi_map_activate(GtkMenuItem *m){
  // menubar item entry was clicked
  GtkWidget *dialog;
  FILE *fpointer;
  char filenameMap[512];
  char singleLine[256];
  char temp[32];
  float DistCenterCloud,rPlus,rMinus;
  int cloudArrayPositioning=0;

  // open map file
  dialog = gtk_file_chooser_dialog_new("Open Map File",GTK_WINDOW(window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_OK,GTK_RESPONSE_OK,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,NULL);
  gtk_widget_show_all(dialog);
  gint resp=gtk_dialog_run(GTK_DIALOG(dialog));
  if(resp==GTK_RESPONSE_OK){
    strcpy(filenameMap,gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog)));
  }else{
    // nothing here
  }
  gtk_widget_destroy(dialog);

  // read file and filter out the data to store in corresponding variables
  fpointer = fopen(filenameMap,"r");
  int counter=0;
  while(!feof(fpointer)){
    gtk_tree_store_append(storeMap,&mapIter,NULL);
    fgets(singleLine, 256, fpointer);
    char *temp1[64];
    char *temp2[64];
    char *temp3[64];
    char *temp4[64];
    fscanf(fpointer,"%s %s %s %s",&temp1,&temp2,&temp3,&temp4);
    mapData[counter][0] = strtof((const char *restrict) temp1,NULL);
    mapData[counter][1] = strtof((const char *restrict) temp2,NULL);
    mapData[counter][2] = strtof((const char *restrict) temp3,NULL);
    mapData[counter][3] = strtof((const char *restrict) temp4,NULL);

    // g_print("%f %f %f %f\n",mapData[counter][0],mapData[counter][1],mapData[counter][2],mapData[counter][3]); // debug

    gtk_tree_store_set(storeMap, &mapIter,0,   (gfloat) mapData[counter][0],-1);
    gtk_tree_store_set(storeMap, &mapIter,1,   (gfloat) mapData[counter][2],-1);
    DistCenterCloud=(220*8.5*sin(mapData[counter][0]*M_PI/180))/(220*sin(mapData[counter][0]*M_PI/180)+mapData[counter][2]);
    gtk_tree_store_set(storeMap, &mapIter,2,   (gfloat) DistCenterCloud,-1);
    rPlus=8.5*cos(mapData[counter][0]*M_PI/180)+sqrt(pow((8.5*cos(mapData[counter][0]*M_PI/180)),2)-pow(8.5,2)+pow(DistCenterCloud,2));
    gtk_tree_store_set(storeMap, &mapIter,3,   (gfloat) rPlus,-1);
    rMinus=8.5*cos(mapData[counter][0]*M_PI/180)-sqrt(pow((8.5*cos(mapData[counter][0]*M_PI/180)),2)-pow(8.5,2)+pow(DistCenterCloud,2));
    gtk_tree_store_set(storeMap, &mapIter,4,   (gfloat) rMinus,-1);
    if(rMinus<=0&&rPlus>=0){ // find a way to display more data than with cancellation of rMinus and rPlus
      // if(1==1){
      gtk_tree_store_set(storeMap, &mapIter,5, (gfloat) rPlus,-1);
      cloudPositionX[cloudArrayPositioning]=rPlus*cos((mapData[counter][0]-90)*M_PI/180);
      gtk_tree_store_set(storeMap, &mapIter,6, (gfloat) cloudPositionX[cloudArrayPositioning], -1);
      cloudPositionY[cloudArrayPositioning]=rPlus*sin((mapData[counter][0]-90)*M_PI/180);
      gtk_tree_store_set(storeMap, &mapIter,7, (gfloat) cloudPositionY[cloudArrayPositioning], -1);
      colorOfMap[cloudArrayPositioning]=(float)mapData[counter][3];
      cloudArrayPositioning++;
    }else{
      gtk_tree_store_set(storeMap, &mapIter,5, (gfloat) 0, -1);
      gtk_tree_store_set(storeMap, &mapIter,6, (gfloat) 0, -1);
      gtk_tree_store_set(storeMap, &mapIter,7, (gfloat) 0, -1);
    }

    counter++;
  }
  fclose(fpointer);

  gtk_widget_queue_draw(mapDrawingArea);
  gtk_widget_show(mapWindow);
}

void on_mi_map_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_map_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> help version
void on_mi_version_activate(GtkMenuItem *m){
  // menubar item entry was clicked
}

void on_mi_version_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_version_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

// Menubar Item -> help help
void on_mi_help_activate(GtkMenuItem *m){
  // menubar item entry was clicked
}

void on_mi_help_select(GtkMenuItem *m){
  // menubar item entry hover on enter
}

void on_mi_help_deselect(GtkMenuItem *m){
  // menubar item entry hover on exit
}

void on_btn_fileAdd_clicked(GtkButton *b){
  // opens the file open dialog
  on_mi_open_activate(NULL);
}

void on_btn_fileDel_clicked(GtkButton *b){
  // deletes selected file from list and clear active preview
  model = gtk_tree_view_get_model(treeViewFiles);
  gtk_tree_store_remove(GTK_TREE_STORE(model),&iter);
  gtk_tree_selection_unselect_all(fileSelection);

  for(int i=0;i<512;i++){
    graphX[i]=0;
    graphY[i]=0;
  }
  gtk_widget_queue_draw(graphDraw);  
  gtk_text_buffer_set_text(textBuffer,(const gchar *) "",(gint) -1);
  gtk_label_set_text((GtkLabel *)lbl_header,(const gchar *)"        --- NO DATA LOADED ---");
  fileListLength--;
}

void on_btn_filePreview_clicked(GtkButton *b){
  // writes data from file to the textbox and draws graph
  reduced = FALSE;
  previewActive=TRUE;
  gtk_widget_set_sensitive( btn_scaleIncrease, TRUE );
  gtk_widget_set_sensitive( btn_scaleDecrease, TRUE );
  // clear out all previous edits
  reduced=FALSE;
  gaussFit=FALSE;
  gtk_widget_set_sensitive( btn_baseline,          TRUE  );
  gtk_widget_set_sensitive( btn_substractBaseline, FALSE );
  gtk_widget_set_sensitive( btn_gaussFit,          FALSE );
  gtk_widget_set_sensitive( btn_saveGauss,         FALSE );
  gtk_widget_set_sensitive( btn_clearCurve,        FALSE );
  for(int i=0;i<4;i++){
    baseX[i]=baseY[i]=0;
  }
  for(int i=0;i<255;i++){
    graphY[i] = graphY[i]+baselineCoordsCopy[i];
    baselineCoords[i]=0;
    baselineCoordsCopy[i]=0;
  }
  
  memset(strLblHeader,0,strlen(strLblHeader));
  gtk_label_set_text((GtkLabel *)lbl_header,(const gchar *)"        --- NO DATA LOADED ---");  
  gchar *fn;
  model = gtk_tree_view_get_model(treeViewFiles);
  gtk_tree_model_get(model,&iter,0,&fn,-1);
  loadTextFile((char **)fn);
}

void on_btnNextSelection_clicked(GtkButton *b){
  // select next item in list
}

void on_btnPrevSelection_clicked(GtkButton *b){
  // select previous item in list
}

void on_fileSelection_changed(GtkWidget *c){
  gchar *fn;
  // check if file is selected or not. if so activate or deactivate filemanagement control buttons
  if(gtk_tree_selection_get_selected(GTK_TREE_SELECTION(c),&model,&iter) == FALSE){
    gtk_widget_set_sensitive( btn_filePreview, FALSE );
    gtk_widget_set_sensitive( btn_fileDel,     FALSE );
    return;
  }else{
    gtk_widget_set_sensitive( btn_filePreview, TRUE );
    gtk_widget_set_sensitive( btn_fileDel,     TRUE );
  }
}

void on_btnTempIncrease_clicked(GtkButton *b){
  // increase temperature scale on graph (y axis)
  scaleModifier++;
  gtk_widget_queue_draw(graphDraw);
  
  // disable button if scaleModifier is 5
  if(scaleModifier<5&&previewActive){
    gtk_widget_set_sensitive( btn_scaleIncrease, TRUE  );
    gtk_widget_set_sensitive( btn_scaleDecrease, TRUE  );
  }else{
    gtk_widget_set_sensitive( btn_scaleIncrease, FALSE );
    gtk_widget_set_sensitive( btn_scaleDecrease, TRUE  );
  }

  // scale labels according to scaleModifier
  if(scaleModifier==1){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "40"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "60"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "80"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "100"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "120"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "140"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "160"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "180"  );
  }else if(scaleModifier==2){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "10"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "30"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "40"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "50"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "60"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "70"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "80"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "90"   );
  }else if(scaleModifier==3){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "6.7"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "13.4" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "20.1" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "26.8" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "33.5" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "40.2" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "46.9" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "53.6" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "60.3" );
  }else if(scaleModifier==4){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "5"    );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "10"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "15"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "25"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "30"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "35"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "40"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "45"   );
  }else if(scaleModifier==5){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "4"    );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "8"    );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "12"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "16"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "24"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "28"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "32"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "36"   );
  }
}

void on_btnTempDecrease_clicked(GtkButton *b){
  // decrease temperature scale on graph (y axis)
  scaleModifier--;
  gtk_widget_queue_draw(graphDraw);

  // disable button if scaleModifier is 1
  if(scaleModifier>1&&previewActive){
    gtk_widget_set_sensitive( btn_scaleDecrease, TRUE  );
    gtk_widget_set_sensitive( btn_scaleIncrease, TRUE  );
  }else{
    gtk_widget_set_sensitive( btn_scaleDecrease, FALSE );
    gtk_widget_set_sensitive( btn_scaleIncrease, TRUE  );
  }

  // scale labels according to scaleModifier
  if(scaleModifier==1){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "40"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "60"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "80"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "100"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "120"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "140"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "160"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "180"  );
  }else if(scaleModifier==2){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "10"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "30"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "40"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "50"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "60"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "70"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "80"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "90"   );
  }else if(scaleModifier==3){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "6.7"  );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "13.4" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "20.1" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "26.8" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "33.5" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "40.2" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "46.9" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "53.6" );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "60.3" );
  }else if(scaleModifier==4){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "5"    );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "10"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "15"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "25"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "30"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "35"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "40"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "45"   );
  }else if(scaleModifier==5){
    gtk_label_set_text( (GtkLabel *) lbl_Temp1, (const gchar *) "4"    );
    gtk_label_set_text( (GtkLabel *) lbl_Temp2, (const gchar *) "8"    );
    gtk_label_set_text( (GtkLabel *) lbl_Temp3, (const gchar *) "12"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp4, (const gchar *) "16"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp5, (const gchar *) "20"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp6, (const gchar *) "24"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp7, (const gchar *) "28"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp8, (const gchar *) "32"   );
    gtk_label_set_text( (GtkLabel *) lbl_Temp9, (const gchar *) "36"   );
  }
}
// print matrix for debugging
int N,M=4;
void printMatrix(const double mat[4][4]){
  for(size_t i=0;i<4;++i){
    for(size_t j=0;j<4;++j){
      printf("%0.15f\t",mat[i][j]);
    }
    printf("\n");
  }
}

// swapping line of matrix
gboolean swapLine(double mat[N][M], unsigned int line1, unsigned int line2){
  if(line1>=N || line2 >=N)
    return FALSE;

  for(size_t i=0;i<M;i++){
    double t = mat[line1][i];
    mat[line1][i]=mat[line2][i];
    mat[line2][i]=t;
  }
  return TRUE;
}

// invert matrix
int N=4;
gboolean invertMatrix(const double mat[N][N], double inv[N][N]){
  double A[N][2*N];
  for(size_t i=0;i<N;++i){
    for(size_t j=0;j<N;++j)
      A[i][j]=mat[i][j];
    for(size_t j=N;j<2*N;++j)
      A[i][j]=(i==j-N)?1.0:0.0;
  }

  for(size_t k=0;k<N-1;++k){
    if(A[k][k]==0.0){
      for(size_t i=k+1;i<N;++i){
        if(A[i][k]!=0.0){
          N=N;
          M=2*N;
          swapLine(A,k,i);
          break;
        }else if(i==N-1){
          return FALSE;
        }
      }
    }
    for(size_t i=k+1;i<N;++i){
      double p = A[i][k]/A[k][k];
      for(size_t j=k;j<2*N;++j){
        A[i][j]-=A[k][j]*p;
      }
    }
  }
  double det = 1.0;
  for(size_t k = 0;k<N;++k)
    det *= A[k][k];
  if(det == 0.0)
    return FALSE;

  for(size_t k=N-1;k>0;--k){
    for(int i = k-1;i>=0;--i){
      double p = A[i][k]/A[k][k];
      for(size_t j = k;j<2*N;++j)
        A[i][j] -= A[k][j]*p;
    }
  }
  for(size_t i=0;i<N;++i){
    const double f = A[i][i];
    for(size_t j=N;j<2*N;++j)
      inv[i][j-N]=A[i][j]/f;
  }
  return TRUE;
}

// multiply Matrices
gboolean multiplyMatrix(double invA[4][4], double B[4],double res[4]){
  double t;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      t+=invA[i][j]*B[j];
    }
    res[i]=t;
    t=0;
  }
  return TRUE;
}

void on_btn_baseline_clicked(GtkButton *b){
  // calculate baseline for spektrum from the previous taken four points
  if(baseX[0]!=0&&baseY[0]!=0&&baseX[1]!=0&&baseY[1]!=0&&baseX[2]!=0&&baseY[2]!=0&&baseX[3]!=0&&baseY[3]!=0){
    for(int i=0;i<4;i++){
      XMatrix[i][0] = pow(baseX[i],3);
      XMatrix[i][1] = pow(baseX[i],2);
      XMatrix[i][2] = baseX[i];
      XMatrix[i][3] = 1;
      YMatrix[i] = baseY[i];
    }
    invertMatrix(XMatrix,XMatrixInv);
    multiplyMatrix(XMatrixInv,YMatrix,baselineMatrix);

    for(int i=0;i<255;i++){
      baselineCoords[i]=baselineMatrix[0]*pow(graphX[i],3)+baselineMatrix[1]*pow(graphX[i],2)+baselineMatrix[2]*graphX[i]+baselineMatrix[3];
    }
    // enable and disable appropriate tool buttons
    gtk_widget_set_sensitive( btn_baseline,          FALSE );
    gtk_widget_set_sensitive( btn_substractBaseline, TRUE  );
    gtk_widget_set_sensitive( btn_gaussFit,          FALSE );
    gtk_widget_set_sensitive( btn_saveGauss,         FALSE );
    gtk_widget_set_sensitive( btn_clearCurve,        TRUE  );
    gtk_widget_queue_draw(graphDraw);
  }else{
    g_print("too few points selected for baseline curve calculation\n");
  }
}

void on_btn_substractBaseline_clicked(GtkButton *b){
  // substract baseline from spectrum
  for(int i=0;i<255;i++){
    graphY[i]-=baselineCoords[i];
    baselineCoordsCopy[i]=baselineCoords[i];
    baselineCoords[i]=0;
  }
  reduced=TRUE;
  gaussFit=TRUE;
  for(int i=0;i<4;i++){
    baseX[i]=baseY[i]=0;
  }
  baseX[0]=-300;
  baseX[1]=300;
  // enable and disable appropriate tool buttons
  gtk_widget_set_sensitive( btn_baseline,          FALSE );
  gtk_widget_set_sensitive( btn_substractBaseline, FALSE );
  gtk_widget_set_sensitive( btn_gaussFit,          TRUE  );
  gtk_widget_set_sensitive( btn_saveGauss,         FALSE );
  gtk_widget_set_sensitive( btn_clearCurve,        TRUE  );
  gtk_widget_queue_draw(graphDraw);
}

double peak=0;
double velocity;

void on_btn_gaussFit_clicked(GtkButton *b){
  double mean=0;
  int dataCounter=0;
  double variance=0;
  double stdDev=0;
  double intensity=0;
  peak=0;

  for(int i=0;i<255;i++){
    gaussBellCurve[i]=0;
  }
  // calculate mean
  for(int i=0;i<255;i++){
    if(graphX[i]>=baseX[0]&&graphX[i]<=baseX[1]){
      if(graphY[i]>=peak){
        peak=graphY[i];
        velocity=graphX[i];
        mean=graphX[i];
      }
    }
  }
  // mean/=dataCounter;

  // calculate variance and therefore standard deviation
  for(int i=0;i<255;i++){
    if(graphX[i]>=baseX[0]&&graphX[i]<=baseX[1]){
      variance+=pow(graphY[i]-mean,2);
      dataCounter++;
    }
  }
  variance/=dataCounter;

  double x1,x2;
  x1=baseX[0];
  x2=baseX[1];
  stdDev=(x2-x1)/3;
  if(stdDev<0)
    stdDev*=-1;
  intensity=peak*sqrt(2*M_PI*pow(stdDev,2));
  // calculate gauss fit
  for(int i=0;i<255;i++){
    gaussBellCurve[i]=(1/(stdDev*sqrt(2*M_PI)))*exp(((-0.5)*pow(((graphX[i]-mean)/(stdDev)),2)));
    gaussBellCurve[i]*=intensity;
  }
  baseX[0]=-300;
  baseX[1]=300;

  gtk_widget_queue_draw(graphDraw);

  // enable and disable appropriate tool buttons
  gtk_widget_set_sensitive( btn_baseline,          FALSE );
  gtk_widget_set_sensitive( btn_substractBaseline, FALSE );
  gtk_widget_set_sensitive( btn_gaussFit,          TRUE  );
  gtk_widget_set_sensitive( btn_saveGauss,         TRUE  );
  gtk_widget_set_sensitive( btn_clearCurve,        TRUE  );
}

double outputData[2000][4];
int outputDataCounter=0;
void on_btn_saveGauss_clicked(GtkButton *b){
  // glon ausgabe
  // g_print("%f\t%f\t%f\t%f\n",gLon,gLat,velocity,peak);
  outputData[outputDataCounter][0]=gLon;
  outputData[outputDataCounter][1]=gLat;
  outputData[outputDataCounter][2]=velocity;
  outputData[outputDataCounter][3]=peak;
  outputDataCounter++;
  // add data here
  gtk_tree_store_clear(treeStoreGauss);
  for(int i=0;i<2000;i++){
    if(outputData[i][0]!=0&&outputData[i][1]!=0&&outputData[i][2]!=0&&outputData[i][3]!=0){
      gtk_tree_store_append(treeStoreGauss,&gaussIter,NULL);
      gtk_tree_store_set(treeStoreGauss,&gaussIter,0,(gfloat) outputData[i][0],-1);
      gtk_tree_store_set(treeStoreGauss,&gaussIter,1,(gfloat) outputData[i][1],-1);
      gtk_tree_store_set(treeStoreGauss,&gaussIter,2,(gfloat) outputData[i][2],-1);
      gtk_tree_store_set(treeStoreGauss,&gaussIter,3,(gfloat) outputData[i][3],-1);
    }
  }
}

void on_btn_clearCurve_clicked(GtkButton *b){
  // clear all baseline and gauss line calculations out
  reduced=FALSE;
  gaussFit=FALSE;
  gtk_widget_set_sensitive( btn_baseline,          TRUE  );
  gtk_widget_set_sensitive( btn_substractBaseline, FALSE );
  gtk_widget_set_sensitive( btn_gaussFit,          FALSE );
  gtk_widget_set_sensitive( btn_saveGauss,         FALSE );
  gtk_widget_set_sensitive( btn_clearCurve,        FALSE );
  for(int i=0;i<4;i++){
    baseX[i]=baseY[i]=0;
  }
  for(int i=0;i<255;i++){
    graphY[i] = graphY[i]+baselineCoordsCopy[i];
    baselineCoords[i]=0;
    baselineCoordsCopy[i]=0;
    gaussBellCurve[i]=0;
  }
  gtk_widget_queue_draw(graphDraw);
}

gboolean on_graphDraw_draw(GtkWidget *widget, cairo_t *cr, gpointer data){
  width = gtk_widget_get_allocated_width(widget);
  height = gtk_widget_get_allocated_height(widget);

  // draw coordinate system
  cairo_set_line_width(cr,1.0);
  cairo_set_source_rgb(cr,0.3,0.3,0.3);
  cairo_move_to(cr, (double) 30.0,                 (double) 0.0           );
  cairo_line_to(cr, (double) 30.0,                 (double) height-30     );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 30.0,                 (double) height-30     );
  cairo_line_to(cr, (double) width,                (double) height-30     );
  cairo_stroke(cr);
  cairo_set_source_rgb(cr,0.75,0.75,0.75);
  cairo_move_to(cr, (double) -300*0.88+width/2+16, (double) 0             );
  cairo_line_to(cr, (double) -300*0.88+width/2+16, (double) height        );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) -200*0.88+width/2+16, (double) 0             );
  cairo_line_to(cr, (double) -200*0.88+width/2+16, (double) height        );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) -100*0.88+width/2+16, (double) 0             );
  cairo_line_to(cr, (double) -100*0.88+width/2+16, (double) height        );
  cairo_stroke(cr);
  cairo_set_source_rgb(cr,0.5,0.5,0.5);
  cairo_move_to(cr, (double) 0*0.88+width/2+16,    (double) 0             );
  cairo_line_to(cr, (double) 0*0.88+width/2+16,    (double) height        );
  cairo_stroke(cr);
  cairo_set_source_rgb(cr,0.7,0.7,0.7);
  cairo_move_to(cr, (double) 100*0.88+width/2+16,  (double) 0             );
  cairo_line_to(cr, (double) 100*0.88+width/2+16,  (double) height        );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 200*0.88+width/2+16,  (double) 0             );
  cairo_line_to(cr, (double) 200*0.88+width/2+16,  (double) height        );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 300*0.88+width/2+16,  (double) 0             );
  cairo_line_to(cr, (double) 300*0.88+width/2+16,  (double) height        );
  cairo_stroke(cr);

  cairo_move_to(cr, (double) 0,                    (double) height-5*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-5*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-10*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-10*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-15*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-15*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-20*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-20*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-25*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-25*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-30*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-30*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-35*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-35*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-40*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-40*7-30 );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,                    (double) height-45*7-30 );
  cairo_line_to(cr, (double) width,                (double) height-45*7-30 );
  cairo_stroke(cr);

  // draw graph
  if(!reduced){
    cairo_set_source_rgb(cr,1.0,0.0,0.0);
  }else{
    cairo_set_source_rgb(cr,0.0,0.5,0.8);
  }
  for(int i=0;i<255;i++){
    if(graphX[i]!=0.0&&graphY[i]!=0.0){
      cairo_move_to(cr, (double) graphX[i]*0.88+width/2+14,   (double) height-graphY[i]*(1.75*scaleModifier)-30   );
      cairo_line_to(cr, (double) graphX[i+1]*0.88+width/2+14, (double) height-graphY[i+1]*(1.75*scaleModifier)-30 );
      cairo_stroke(cr);
    }else{
      return FALSE;
    }
  }

  // draw area for baseline reduction
  if(!reduced){
    static const double dashedLine[] = {1.0};
    for(int i=0;i<4;i++){
      if(baseX[i]!=0&&baseY[i]!=0){
        cairo_set_source_rgb(cr,0.7,0.0,1.0);
        cairo_arc(cr, (double) baseX[i]*0.88+width/2+14, (double) height-baseY[i]*(1.75*scaleModifier)-30, 3, 0, 2 * M_PI );
        cairo_fill(cr);
        cairo_stroke(cr);
        cairo_set_source_rgb(cr,1.0,0.0,1.0);
        cairo_set_dash(cr,dashedLine,1,0);
        cairo_move_to(cr, (double) baseX[i]*0.88+width/2+14, (double) 0      );
        cairo_line_to(cr, (double) baseX[i]*0.88+width/2+14, (double) height );
        cairo_stroke(cr);
      }
    }
  }
  cairo_set_dash(cr,0,0,0);

  // draw baseline
  cairo_set_source_rgb(cr,0.0,0.0,1.0);
  for(int i=0;i<255;i++){
    cairo_move_to(cr, (double) graphX[i]*0.88+width/2+14,   (double) height-baselineCoords[i]*(1.75*scaleModifier)-30   );
    cairo_line_to(cr, (double) graphX[i+1]*0.88+width/2+14, (double) height-baselineCoords[i+1]*(1.75*scaleModifier)-30 );
    cairo_stroke(cr);
  }

  // draw area for gauss fit detection
  if(gaussFit){
    static const double dashedLine[] = {1.0};
    for(int i=0;i<2;i++){
      cairo_set_source_rgb(cr,0.0,0.5,0.0);
      cairo_set_dash(cr,dashedLine,1,0);
      cairo_move_to(cr, (double) baseX[i]*0.88+width/2+14, (double) 0      );
      cairo_line_to(cr, (double) baseX[i]*0.88+width/2+14, (double) height );
      cairo_stroke(cr);
    }
  }
  cairo_set_dash(cr,0,0,0);

  // draw gauss bell curve
  if(gaussFit){
    cairo_set_line_width(cr,2.0);
    cairo_set_source_rgb(cr,0.3,1.0,0.3);
    for(int i=0;i<255;i++){
      cairo_move_to(cr, (double) graphX[i]*0.88+width/2+14,   (double) height-gaussBellCurve[i]*(1.75*scaleModifier)-30   );
      cairo_line_to(cr, (double) graphX[i+1]*0.88+width/2+14, (double) height-gaussBellCurve[i+1]*(1.75*scaleModifier)-30 );
      cairo_stroke(cr);
    }
  }
  
  return FALSE;
}

void on_btn_closeMap_clicked(GtkButton *b){
  gtk_widget_hide(mapWindow);
}

void on_btn_export3dMap_clicked(GtkButton *b){

}

void on_btn_exportMapImage_clicked(GtkButton *b){

}

// maybe make the following interactive for zooming in and out 
int colorModifier = 80;
int zoomModifier = 10;

gboolean on_mapDrawingArea_draw(GtkWidget *widget, cairo_t *cr, gpointer data){
  width = gtk_widget_get_allocated_width(widget);
  height = gtk_widget_get_allocated_height(widget);
  
  // draw coordinate system
  cairo_set_line_width(cr,1.0);
  cairo_set_source_rgb(cr,0.0,0.0,0.0);
  cairo_move_to(cr, (double) width/2, (double) 0        );
  cairo_line_to(cr, (double) width/2, (double) height   );
  cairo_stroke(cr);
  cairo_move_to(cr, (double) 0,       (double) height/2 );
  cairo_line_to(cr, (double) width,   (double) height/2 );
  cairo_stroke(cr);

  // draw hi clouds
  int temperatureIter=0;
  for(int i=0;i<2000;i++){
    if(i%256==0){
      temperatureIter++;
    }
    // find out why temperature is displayed crappy
    if(cloudPositionX[i]==0 && cloudPositionY[i]==0)
      i=1999;

    cairo_set_source_rgb(cr,0.3,0.6,colorOfMap[i]/100);
    cairo_arc(cr, (double) width/2+cloudPositionX[i]*zoomModifier, (double) height/2+cloudPositionY[i]*zoomModifier , 3, 0, 2 * M_PI );
    cairo_fill(cr);
    cairo_stroke(cr);

  }
  return FALSE;
}

void on_mapSelection_changed(GtkWidget *c){

}

gboolean on_graphDraw_button_press_event(GtkWidget *widget, GdkEventButton *event, gpointer user_data){
  width = gtk_widget_get_allocated_width(widget);
  height = gtk_widget_get_allocated_height(widget);

  if(!gaussFit&&  gtk_widget_get_sensitive(btn_baseline)){ // get position of points for baseline curve
    reduced=FALSE;
    baseX[clickCounter]=(event->x-width/2-16)*1.136;
    baseY[clickCounter]=((event->y-height)*(-1)-30)*0.573689/scaleModifier;

    if(clickCounter==3){
      clickCounter=0;
    }else{
      clickCounter++;
    }
    gtk_widget_queue_draw(widget);
  }else if(gaussFit){ // get position of points for gauss fit
    baseX[clickCounter]=(event->x-width/2-16)*1.136;
    baseY[clickCounter]=((event->y-height)*(-1)-30)*0.573689/scaleModifier;
    if(clickCounter==1){
      clickCounter=0;
    }else{
      clickCounter++;
    }
    gtk_widget_queue_draw(widget);
  }
  return FALSE;
}

void on_gaussSelection_changed(GtkWidget *c){

}

void on_btn_deleteGaussCoords_clicked(GtkButton *b){

}

void on_btn_deleteAllGaussCoords_clicked(GtkButton *b){
  for(int i=0;i<2000;i++){
    outputData[i][0]=0;
    outputData[i][1]=0;
    outputData[i][2]=0;
    outputData[i][3]=0;
    outputDataCounter=0;
  }
  gtk_tree_store_clear(treeStoreGauss);
}

void on_btn_exportGaussCoords_clicked(GtkButton *b){
  GtkWidget *dialog;
  char filenameSave[512];
  FILE *fpointer;

  dialog = gtk_file_chooser_dialog_new("Expxort Data",GTK_WINDOW(window),GTK_FILE_CHOOSER_ACTION_SAVE,GTK_STOCK_OK,GTK_RESPONSE_OK,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,NULL);
  gtk_widget_show_all(dialog);
  gint resp=gtk_dialog_run(GTK_DIALOG(dialog));
  if(resp==GTK_RESPONSE_OK){
    strcpy(filenameSave,(const char *restrict)gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog)));
    g_print("%s\n",filenameSave); // debug
  }else{
    // nothing here
  }
  gtk_widget_destroy(dialog);

  fpointer = fopen(filenameSave,"w");
  for(int i=0;i<2000;i++){
    if(outputData[i][0]!=0&&outputData[i][1]!=0&&outputData[i][2]!=0&&outputData[i][3]!=0){
      fprintf(fpointer,"%f %f %f %f\n",outputData[i][0],outputData[i][1],outputData[i][2],outputData[i][3]);
    }
  }
  fclose(fpointer);
}
