// aquat.h -- Plotting class using AquaTerm
// language: C++
// author  : Paris Smaragdis

#ifndef __AQUAT_H__
#define __AQUAT_H__

#ifdef AQUA

#include <cmath>
#include <algorithm>

extern "C"{
	#include "aquaterm.h"
	int _aq_w = 400, _aq_h = 200;
}

// Image colormaps
void set_color( float c, int m, float p)
{
	float r = 0, g = 0, b = 0;
	switch( m){
		case 1:
			// Get RGB values
			r = 1;
			if( c < .36)
				r = .004 + (c/.36) * (1.-.004);
			g = 0;
			if( c > .74)
				g = 1;
			else if( c > .38)
				g = (c-.38)/(.74-.38);
			b = 0;
			if( c > .74)
				b = (c-.74)/(1-.74);
		break;
		default:
			r = g = b = c;
	}

	// Clip to 0,1 interval
	r = std::min( 1.f, std::max( 0.f, r));
	g = std::min( 1.f, std::max( 0.f, g));
	b = std::min( 1.f, std::max( 0.f, b));

	// Set the color
	aqtSetColor( std::pow( r, p), std::pow( g, p), std::pow( b, p));
}

// Open an aquaterm window by id
void aq_window( int id = 1, const char *c = NULL)
{
	// Initialize if we haven't done so
	static bool init = false;
	if( !init){
		aqtInit();
		init = true;
	}

	// Open the window with the provided id
	aqtOpenPlot( id);
	if( c)
		aqtSetPlotTitle( c);
	aqtSetPlotSize( _aq_w, _aq_h);
}

// Clear the current window
void aq_clear()
{
	aqtClearPlot();
}

// Plot a sequence with a line
template <class T>
void aq_plot( const T *x, int l, int r = 0, int g = 0, int b = 0)
{
	aqtClearPlot();
	
	// Find the extrema
	T mn = x[0], mx = x[0];
	for( int i = 0 ; i < l ; i++){
		mn = min( mn, x[i]);
		mx = max( mx, x[i]);
	}
	
	// Draw the line
	aqtSetLinestyleSolid();
	aqtSetColor( r, g, b);
	T ln = .9*_aq_h/(mx - mn);
	for( int i = 0 ; i < l-1 ; i++){
		aqtMoveTo( _aq_w*float(i)/(l-1), .05*_aq_h+ln*(x[i]-mn));
		aqtAddLineTo( _aq_w*float(i+1)/(l-1), .05*_aq_h+ln*(x[i+1]-mn));
	}
	aqtRenderPlot();
}

// Plot multiple sequences with different colored lines
template <class T>
void aq_mplot( const T *x, int n, int l)
{
	const float cr[] = {1,0,0,1,1,0,0};
	const float cb[] = {0,1,0,1,0,1,0};
	const float cg[] = {0,0,1,0,1,1,0};

	aqtClearPlot();
	
	// Find the extrema
	T mn = x[0], mx = x[0];
	for( int i = 0 ; i < n*l ; i++){
		mn = min( mn, x[i]);
		mx = max( mx, x[i]);
	}
	
	// Draw the line
	for( int j = 0 ; j < n ; j++){
		aqtSetLinestyleSolid();
		aqtSetColor( cr[j%7], cg[j%7], cb[j%7]);
		T ln = .9*_aq_h/(mx - mn);
		for( int i = 0 ; i < l-1 ; i++){
			aqtMoveTo( _aq_w*float(i)/(l-1), .05*_aq_h+ln*(x[j*l+i]-mn));
			aqtAddLineTo( _aq_w*float(i+1)/(l-1), .05*_aq_h+ln*(x[j*l+i+1]-mn));
		}
	}
	aqtRenderPlot();
}

// Plot a point sequence
template <class T>
void aq_plotp( const T *x, int l, int sz = 2, int r = 0, int g = 0, int b = 0)
{
	aqtClearPlot();

	// Find the extrema
	T mn = x[0], mx = x[0];
	for( int i = 0 ; i < l ; i++){
		mn = min( mn, x[i]);
		mx = max( mx, x[i]);
	}

	// Draw the line
	aqtSetLinestyleSolid();
	aqtSetLineCapStyle( AQTRoundLineCapStyle);
	aqtSetColor( r, g, b);
	T ln = .9*_aq_h/(mx - mn);
	for( int i = 0 ; i < l-1 ; i++)
		aqtAddFilledRect( _aq_w*float(i)/(l-1), .05*_aq_h+ln*(x[i]-mn), sz, sz);
	aqtRenderPlot();
}

template <class T>
void aq_image( const T *x, int h, int w, int m = 0, float p = 1)
{
	// Clean up the previous plot
	aqtClearPlot();

	// Find the extrema
	T mn = x[0], mx = x[0];
	for( int i = 0 ; i < w*h ; i++){
		mn = min( mn, x[i]);
		mx = max( mx, x[i]);
	}

	// Make the plot
	for( int i = 0 ; i < h ; i++)
		for( int j = 0 ; j < w ; j++){
			float c = (x[i+h*j]-mn) / (mx-mn);
			set_color( c, m, p);
			aqtAddFilledRect( _aq_w*float(j)/(w), _aq_h*float(i)/(h), float(_aq_w)/w, float(_aq_h)/h);
		}
	aqtRenderPlot();
}

// Add some text
void aq_text( char *s, double x, double y)
{
	aqtSetFontname( "Monaco");
	aqtSetFontsize( 10);
	aqtAddLabel( s, _aq_w*x, _aq_h*y, 0, AQTAlignCenter);
	aqtRenderPlot();
}

// Close a window
void aq_close()
{
//	aqtClosePlot();
//	aqtTerminate();
}

#else

void aq_window( int id = 1, const char *c = NULL) {}
void aq_clear() {}
template <class T> void aq_plot( const T *x, int l, int r = 0, int g = 0, int b = 0) {}
template <class T> void aq_plotp( const T *x, int l, int sz = 2, int r = 0, int g = 0, int b = 0) {}
template <class T> void aq_image( const T *x, int h, int w, int m = 0, float p = 1) {}
void aq_text( char *s, double x, double y) {}
void aq_close() {}

#endif


#endif
