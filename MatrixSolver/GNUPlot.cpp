#include <xlocale.h>
#include "GNUPlot.hpp"

#include <iostream>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <unistd.h>

using namespace std;

GNUPlot::GNUPlot(void) : m_pPlot(NULL), m_nPlots(0), m_nTemp(0) 
{ 
    if (getenv("DISPLAY") == NULL)
        std::cerr << "Cannot Find DISPLAY Variable!" << std::endl;
    
    if ( GetProgramPath("gnuplot") == NULL )
        std::cerr << "Cannot Find gnuplot in the PATH" << std::endl;
    
    SetStyle("lines");
    
    m_pPlot = popen("gnuplot -persist", "w");
    
    if (m_pPlot == NULL)
        std::cerr << "Error Starting gnuplot" << std::endl;
}

GNUPlot::GNUPlot(const GNUPlot& plot)
{
    std::cerr << "Calling Copy Constructor for GNUPlot Class. " << std::endl;
    m_pPlot   = plot.m_pPlot;
    strcpy(m_szStyle, plot.m_szStyle);
    m_nPlots  = plot.m_nPlots;
    m_nTemp   = plot.m_nTemp;
    for (unsigned int i = 0; i < GNUPLOT_MAX_TEMP_FILES; ++i)
        strcpy(m_szTemp[i], plot.m_szTemp[i]);
}

GNUPlot::~GNUPlot(void)
{
    Close();
}

void GNUPlot::Close(void)
{
    unsigned int i;
    
    if ( pclose(m_pPlot) == -1 )
    {	
        std::cerr << "Problem Closing gnuplot" << std::endl;
        return;
    }
    
    if (m_nTemp)
    {
        for (i = 0; i < m_nTemp; ++i)
            remove(m_szTemp[i]);
    }
}

void GNUPlot::Command(const char* szCmd, ...)
{
    va_list vaList;
    char szLocal[GNUPLOT_COMMAND_SIZE];
    
    va_start(vaList, szCmd);
    vsprintf(szLocal, szCmd, vaList);
    va_end(vaList);
    
    strcat(szLocal, "\n");
    
    fprintf(m_pPlot, "%s", szLocal);
    fflush(m_pPlot);
}

void GNUPlot::SetStyle(const char* szStyle)
{
    if ( strcmp(szStyle, "lines"        ) && strcmp(szStyle, "points"   ) &&
        strcmp(szStyle, "linespoints"  ) && strcmp(szStyle, "impulses" ) &&
        strcmp(szStyle, "dots"         ) && strcmp(szStyle, "steps"    ) &&
        strcmp(szStyle, "errorbars"    ) && strcmp(szStyle, "boxes"    ) &&
        strcmp(szStyle, "boxerrorbars" ) )
    {
        std::cerr << "Warning (GNUPlot): Unknown Style Requested, Using 'points' " << std::endl;
        strcpy(m_szStyle, "points");
    }
    else
    {
        strcpy(m_szStyle, szStyle);
    }
}

void GNUPlot::SetLabelX(const char* szLabel)
{
    char szCmd[GNUPLOT_COMMAND_SIZE];
    
    sprintf(szCmd, "set xlabel '%s'", szLabel);
    Command(szCmd);
}

void GNUPlot::SetLabelY(const char* szLabel)
{
    char szCmd[GNUPLOT_COMMAND_SIZE];
    
    sprintf(szCmd, "set ylabel '%s'", szLabel);
    Command(szCmd);
}

void GNUPlot::ResetPlot(void)
{
    unsigned int i;
    if (m_nTemp)
    {
        for (i = 0; i < m_nTemp; ++i)
            remove(m_szTemp[i]);
    }
    
    m_nTemp  = 0;
    m_nPlots = 0;
}

void GNUPlot::Plot(double *pX, unsigned int nSize, const char* szTitle)
{
    unsigned int i;
    int fdTemp;
    char szName[128];
    char szCmd[GNUPLOT_COMMAND_SIZE];
    char szLine[GNUPLOT_COMMAND_SIZE];
    
    if ( (m_pPlot == NULL) || (pX == NULL) || (nSize < 1) )
        return;
    
    if (m_nTemp == GNUPLOT_MAX_TEMP_FILES - 1)
    {
        std::cerr << "Maximum # of Temp Files Reached (" << GNUPLOT_MAX_TEMP_FILES << "): Cannot Open More." << std::endl;
        return;
    }
    
    sprintf(szName, "%s/gnuplot-XXXXXX", GNUPLOT_TEMP_DIR);
    if ( ( fdTemp = mkstemp(szName) ) == -1 )
    {
        std::cerr << "Error Creating Temporary File for Plotting. " << std::endl;
        return;
    }
    
    strcpy(m_szTemp[m_nTemp++], szName);
    
    for (i = 0; i < nSize; ++i)
    {
        sprintf(szLine, "%g\n", pX[i]);
        write(fdTemp, szLine, strlen(szLine) );
    }
    close(fdTemp);
    
    if (m_nPlots > 0)
        strcpy(szCmd, "replot");
    else
        strcpy(szCmd, "plot");
    
    if (szTitle == NULL)
        sprintf(szLine, "%s '%s' with %s", szCmd, szName, m_szStyle);
    else
        sprintf(szLine, "%s '%s' title '%s' with %s", szCmd, szName, szTitle, m_szStyle);
    
    Command(szLine);
    m_nPlots++;
}

char* GNUPlot::GetProgramPath(const char* szName)
{
    unsigned int i;
    unsigned int j;
    size_t nLength;
    char* szPath;
    
    static char szBuffer[GNUPLOT_PATH_MAX_NAMES];
    
    sprintf(szBuffer, "./%s", szName);
    if ( access(szBuffer, 00) == 0)
    {
        sprintf(szBuffer, ".");
        return (szBuffer);
    }
    
    szBuffer[0] = 0;
    
    szPath = getenv("PATH");
    
    if (szPath != NULL)
        for (i = 0; szPath[i]; )
        {
            for (j = i; (szPath[j]) && (szPath[j] != ':'); ++j); /* Do Nothing! */
            nLength = j - i;
            strncpy(szBuffer, szPath + i, nLength);
            if (nLength == 0)
                szBuffer[nLength++] = char(46);
            szBuffer[nLength++] = char(47);
            strcpy(szBuffer + nLength, szName);
            if ( access(szBuffer, 00) == 0) /* Found GNUPlot */
                break;
            szBuffer[0] = 0;
            i = j;
            if (szPath[i] == char(58)) 
                ++i;
        }
        else
            std::cerr << "PATH variable not set." << std::endl;
        
        if (szBuffer[0] == 0)
            return (NULL);
        
        nLength = strlen(szBuffer) - 1;
        
        while (szBuffer[nLength] != char(47))
            szBuffer[nLength--] = 0;
        
        szBuffer[nLength] = 0;
        
        return (szBuffer);
}

void GNUPlot::Plot(double* pX, double* pY, unsigned int nSize, const char* szTitle)
{
    unsigned int i;
    int fdTemp;
    char szName[128];
    char szCmd[GNUPLOT_COMMAND_SIZE];
    char szLine[GNUPLOT_COMMAND_SIZE];
    
    if ( (m_pPlot == NULL) || (pX == NULL) || (nSize < 1) )
        return;
    
    if (m_nTemp == GNUPLOT_MAX_TEMP_FILES - 1)
    {
        std::cerr << "Maximum # of Temp Files Reached (" << GNUPLOT_MAX_TEMP_FILES << "): Cannot Open More." << std::endl;
        return;
    }
    
    sprintf(szName, "%s/gnuplot-i-XXXXXX", GNUPLOT_TEMP_DIR);
    if ( ( fdTemp = mkstemp(szName) ) == -1 )
    {
        std::cerr << "Error Creating Temporary File for Plotting. " << std::endl;
        return;
    }
    
    strcpy(m_szTemp[m_nTemp++], szName);
    
    for (i = 0; i < nSize; ++i)
    {
        sprintf(szLine, "%g %g\n", pX[i], pY[i]);
        write(fdTemp, szLine, strlen(szLine) );
    }
    close(fdTemp);
    
    if (m_nPlots > 0)
        strcpy(szCmd, "replot");
    else
        strcpy(szCmd, "plot");
    
    if (szTitle == NULL)
        sprintf(szLine, "%s '%s' with %s", szCmd, szName, m_szStyle);
    else
        sprintf(szLine, "%s '%s' title '%s' with %s", szCmd, szName, szTitle, m_szStyle);
    
    Command(szLine);
    m_nPlots++;
}

void GNUPlot::PlotSlope(double a, double b, const char* szTitle)
{
    char szCmd[GNUPLOT_COMMAND_SIZE];
    
    if (m_nPlots > 0)
        sprintf(szCmd, "replot %g * x + %g title '%s' with %s", a, b, szTitle, m_szStyle);
    else
        sprintf(szCmd, "plot %g * x + %g title '%s' with %s", a, b, szTitle, m_szStyle);
    
    Command(szCmd);
    m_nPlots++;
}

void GNUPlot::PlotEquation(const char* szEquation, const char* szTitle)
{
    char szCmd[GNUPLOT_COMMAND_SIZE];
    
    if (m_nPlots > 0)
        sprintf(szCmd, "replot %s title '%s' with %s", szEquation, szTitle, m_szStyle);
    else
        sprintf(szCmd, "plot %s title '%s' with %s", szEquation, szTitle, m_szStyle);
    
    Command(szCmd);
    m_nPlots++;
}

void GNUPlot::Plot(double* pX, double* pY, double** pZ, unsigned int nSizeX, unsigned int nSizeY, const char* szTitle)
{
    unsigned int i, j;
    int fdTemp;
    char szName[128];
    char szCmd[GNUPLOT_COMMAND_SIZE];
    char szLine[GNUPLOT_COMMAND_SIZE];
    
    if ( (m_pPlot == NULL) || (pX == NULL) || (nSizeX < 1) || (nSizeY < 1) )
        return;
    
    if (m_nTemp == GNUPLOT_MAX_TEMP_FILES - 1)
    {
        std::cerr << "Maximum # of Temp Files Reached (" << GNUPLOT_MAX_TEMP_FILES << "): Cannot Open More." << std::endl;
        return;
    }
    
    sprintf(szName, "%s/gnuplot-i-XXXXXX", GNUPLOT_TEMP_DIR);
    if ( ( fdTemp = mkstemp(szName) ) == -1 )
    {
        std::cerr << "Error Creating Temporary File for Plotting. " << std::endl;
        return;
    }
    
    strcpy(m_szTemp[m_nTemp++], szName);
    
    for (i = 0; i < nSizeX; ++i)
    {
        for (j = 0; j < nSizeY; ++j)
        {
            sprintf(szLine, "%g %g %g\n", pX[i], pY[j], pZ[i][j]);
            write(fdTemp, szLine, strlen(szLine) );
        }
        sprintf(szLine, "\n");
        write(fdTemp, szLine, strlen(szLine) );
    }
    close(fdTemp);
    
    if (m_nPlots > 0)
        strcpy(szCmd, "sreplot");
    else
        strcpy(szCmd, "splot");
    
    if (szTitle == NULL)
        sprintf(szLine, "%s '%s' with %s", szCmd, szName, m_szStyle);
    else
    {
        sprintf(szLine, "set title '%s' ", szTitle);
        Command(szLine);
        sprintf(szLine, "%s '%s' title '%s' with %s", szCmd, szName, szTitle, m_szStyle);
    }
    
    Command(szLine);
    m_nPlots++;
}

void GNUPlot::PlotContour(double* pX, double* pY, double** pZ, unsigned int nSizeX, unsigned int nSizeY, const char* szTitle)
{
    Command("set contour\n");
    Command("set parametric\n");
    Command("set nosurface\n");
    Command("set view 0,0,1\n");
    Command("set cntrparam levels 10\n");
    Command("set grid\n");
    Plot(pX, pY, pZ, nSizeX, nSizeY, szTitle);
}

void GNUPlot::PlotSurface(double* pX, double* pY, double** pZ, unsigned int nSizeX, unsigned int nSizeY, const char* szTitle)
{
    Command("set contour\n");
    Command("set parametric\n");
    Command("set view 1,1,1\n");
    Command("set cntrparam levels 10\n");
    Plot(pX, pY, pZ, nSizeX, nSizeY, szTitle);
}

void GNUPlot::PlotLevelSet(double* pX, double* pY, double** pZ, unsigned int nSizeX, unsigned int nSizeY, const char* szTitle)
{
    Command("set data style lines");
    Command("set contour\n");
    Command("set parametric\n");
    Command("set nosurface\n");
    Command("set view 0,0,1\n");
    Command("set cntrparam levels 1\n");
    Command("set cntrparam levels discrete 0.0\n");
    Command("set grid\n");
    Plot(pX, pY, pZ, nSizeX, nSizeY, szTitle);
}

void GNUPlot::SetXTics(const double Size)
{
    char szLine[GNUPLOT_COMMAND_SIZE];
    sprintf(szLine, "set xtics %g\n", Size);
    Command(szLine);
}

void GNUPlot::SetYTics(const double Size)
{
    char szLine[GNUPLOT_COMMAND_SIZE];
    sprintf(szLine, "set ytics %g\n", Size);
    Command(szLine);
}
