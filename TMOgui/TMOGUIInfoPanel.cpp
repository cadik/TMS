
#include "TMOGUIInfoPanel.h"
#include "TMOGUILineResizer.h"

#include <qlabel.h>
#include <qgroupbox.h>
#include <qtooltip.h>


TMOGUIInfoPanel::TMOGUIInfoPanel( QWidget* parent, const char* name )
    : QWidget( parent, name )
{
    groupBoxStat = new QGroupBox( this, "groupBoxStat" );
    groupBoxStat->setGeometry( QRect( 290, 10, 130, 140 ) );

    pAbove = new QLabel( groupBoxStat, "pAbove" );
    pAbove->setGeometry( QRect( 70, 50, 49, 20 ) );

    textLabelBelow = new QLabel( groupBoxStat, "textLabelBelow" );
    textLabelBelow->setGeometry( QRect( 20, 30, 49, 20 ) );
    QFont textLabelBelow_font(  textLabelBelow->font() );
    textLabelBelow_font.setBold( TRUE );
    textLabelBelow->setFont( textLabelBelow_font ); 

    textLabelAbove = new QLabel( groupBoxStat, "textLabelAbove" );
    textLabelAbove->setGeometry( QRect( 20, 50, 49, 20 ) );
    QFont textLabelAbove_font(  textLabelAbove->font() );
    textLabelAbove_font.setBold( TRUE );
    textLabelAbove->setFont( textLabelAbove_font ); 

    textLabelVisible = new QLabel( groupBoxStat, "textLabelVisible" );
    textLabelVisible->setGeometry( QRect( 20, 70, 49, 20 ) );
    QFont textLabelVisible_font(  textLabelVisible->font() );
    textLabelVisible_font.setBold( TRUE );
    textLabelVisible->setFont( textLabelVisible_font ); 

    textLabelExtreme = new QLabel( groupBoxStat, "textLabelExtreme" );
    textLabelExtreme->setGeometry( QRect( 20, 100, 49, 20 ) );
    QFont textLabelExtreme_font(  textLabelExtreme->font() );
    textLabelExtreme_font.setBold( TRUE );
    textLabelExtreme->setFont( textLabelExtreme_font ); 

    pExtreme = new QLabel( groupBoxStat, "pExtreme" );
    pExtreme->setGeometry( QRect( 70, 100, 49, 20 ) );

    pVisible = new QLabel( groupBoxStat, "pVisible" );
    pVisible->setGeometry( QRect( 70, 70, 49, 20 ) );

    pBelow = new QLabel( groupBoxStat, "pBelow" );
    pBelow->setGeometry( QRect( 70, 30, 49, 20 ) );

    groupBoxHisto = new QGroupBox( this, "groupBoxHisto" );
    groupBoxHisto->setGeometry( QRect( 10, 10, 270, 140 ) );

    pGamma = new QLabel( groupBoxHisto, "pGamma" );
    pGamma->setGeometry( QRect( 80, 70, 180, 20 ) );

    pAverage = new QLabel( groupBoxHisto, "pAverage" );
    pAverage->setGeometry( QRect( 80, 100, 180, 20 ) );

    pWhite = new QLabel( groupBoxHisto, "pWhite" );
    pWhite->setGeometry( QRect( 80, 50, 180, 20 ) );

    pBlack = new QLabel( groupBoxHisto, "pBlack" );	
    pBlack->setGeometry( QRect( 80, 30, 180, 20 ) );

    textLabelBlack = new QLabel( groupBoxHisto, "textLabelBlack" );
    textLabelBlack->setGeometry( QRect( 20, 30, 49, 20 ) );
    QFont textLabelBlack_font(  textLabelBlack->font() );
    textLabelBlack_font.setBold( TRUE );
    textLabelBlack->setFont( textLabelBlack_font ); 

    textLabelWhite = new QLabel( groupBoxHisto, "textLabelWhite" );
    textLabelWhite->setGeometry( QRect( 20, 50, 49, 20 ) );
    QFont textLabelWhite_font(  textLabelWhite->font() );
    textLabelWhite_font.setBold( TRUE );
    textLabelWhite->setFont( textLabelWhite_font ); 

    textLabelGammma = new QLabel( groupBoxHisto, "textLabelGammma" );
    textLabelGammma->setGeometry( QRect( 20, 70, 49, 20 ) );
    QFont textLabelGammma_font(  textLabelGammma->font() );
    textLabelGammma_font.setBold( TRUE );
    textLabelGammma->setFont( textLabelGammma_font ); 

    textLabelAvg = new QLabel( groupBoxHisto, "textLabelAvg" );
    textLabelAvg->setGeometry( QRect( 20, 100, 49, 20 ) );
    QFont textLabelAvg_font(  textLabelAvg->font() );
    textLabelAvg_font.setBold( TRUE );
    textLabelAvg->setFont( textLabelAvg_font ); 

    groupBoxLocalTool = new QGroupBox( this, "groupBoxLocalTool" );
    groupBoxLocalTool->setGeometry( QRect( 430, 10, 510, 140 ) );

    textLabelAvgLum = new QLabel( groupBoxLocalTool, "textLabelAvgLum" );
    textLabelAvgLum->setGeometry( QRect( 20, 50, 95, 20 ) );
    QFont textLabelAvgLum_font(  textLabelAvgLum->font() );
    textLabelAvgLum_font.setBold( TRUE );
    textLabelAvgLum->setFont( textLabelAvgLum_font ); 
	QToolTip::add(textLabelAvgLum, "Average value of luminance in the area");

    pMaxLum = new QLabel( groupBoxLocalTool, "pMaxLum" );
    pMaxLum->setGeometry( QRect( 120, 70, 150, 20 ) );

    pMinLum = new QLabel( groupBoxLocalTool, "pMinLum" );
    pMinLum->setGeometry( QRect( 120, 90, 150, 20 ) );

    textLabelAvgCol = new QLabel( groupBoxLocalTool, "textLabelAvgCol" );
    textLabelAvgCol->setGeometry( QRect( 20, 110, 70, 20 ) );
    QFont textLabelAvgCol_font(  textLabelAvgCol->font() );
    textLabelAvgCol_font.setBold( TRUE );
    textLabelAvgCol->setFont( textLabelAvgCol_font ); 
	QToolTip::add(textLabelAvgCol, "Average value of color in the area <red, green, blue>");

    textLabelMinLum = new QLabel( groupBoxLocalTool, "textLabelMinLum" );
    textLabelMinLum->setGeometry( QRect( 20, 90, 86, 20 ) );
    QFont textLabelMinLum_font(  textLabelMinLum->font() );
    textLabelMinLum_font.setBold( TRUE );
    textLabelMinLum->setFont( textLabelMinLum_font ); 
	QToolTip::add(textLabelMinLum, "Minimum value of luminance in the area");

    textLabelMaxLum = new QLabel( groupBoxLocalTool, "textLabelMaxLum" );
    textLabelMaxLum->setGeometry( QRect( 20, 70, 96, 20 ) );
    QFont textLabelMaxLum_font(  textLabelMaxLum->font() );
    textLabelMaxLum_font.setBold( TRUE );
    textLabelMaxLum->setFont( textLabelMaxLum_font ); 
	QToolTip::add(textLabelMaxLum, "Maximum value of luminance in the area");

    textLabelNeighb = new QLabel( groupBoxLocalTool, "textLabelNeighb" );
    textLabelNeighb->setGeometry( QRect( 20, 20, 89, 20 ) );
    QFont textLabelNeighb_font(  textLabelNeighb->font() );
    textLabelNeighb_font.setBold( TRUE );
    textLabelNeighb_font.setItalic( TRUE );
    textLabelNeighb_font.setUnderline( TRUE );
    textLabelNeighb->setFont( textLabelNeighb_font );
	QToolTip::add(textLabelNeighb, "Data from the place determined by the shape (circle/square)");

    ptextLabelCursor = new QLabel( groupBoxLocalTool, "ptextLabelCursor" );
    ptextLabelCursor->setGeometry( QRect( 280, 20, 89, 20 ) );
    QFont ptextLabelCursor_font(  ptextLabelCursor->font() );
    ptextLabelCursor_font.setBold( TRUE );
    ptextLabelCursor_font.setItalic( TRUE );
    ptextLabelCursor_font.setUnderline( TRUE );
    ptextLabelCursor->setFont( ptextLabelCursor_font ); 
	QToolTip::add(ptextLabelCursor, "Data from the place below the cursor");

    ptextLabelColor = new QLabel( groupBoxLocalTool, "ptextLabelColor" );
    ptextLabelColor->setGeometry( QRect( 280, 70, 40, 20 ) );
    QFont ptextLabelColor_font(  ptextLabelColor->font() );
    ptextLabelColor_font.setBold( TRUE );
    ptextLabelColor->setFont( ptextLabelColor_font ); 
	QToolTip::add(ptextLabelColor, "Color under the cursor");

    pAvgLum = new QLabel( groupBoxLocalTool, "pAvgLum" );
    pAvgLum->setGeometry( QRect( 120, 50, 150, 20 ) );

    ptextLabelLum = new QLabel( groupBoxLocalTool, "ptextLabelLum" );
    ptextLabelLum->setGeometry( QRect( 280, 50, 70, 20 ) );
    QFont ptextLabelLum_font(  ptextLabelLum->font() );
    ptextLabelLum_font.setBold( TRUE );
    ptextLabelLum->setFont( ptextLabelLum_font ); 
	QToolTip::add(ptextLabelLum, "Luminance under the cursor");

    pLum = new QLabel( groupBoxLocalTool, "pLum" );
    pLum->setGeometry( QRect( 350, 50, 150, 20 ) );

    pColor = new QLabel( groupBoxLocalTool, "pColor" );
    pColor->setGeometry( QRect( 350, 70, 150, 20 ) );

    pAvgCol = new QLabel( groupBoxLocalTool, "pAvgCol" );
    pAvgCol->setGeometry( QRect( 120, 110, 150, 20 ) );

	line1 = new TMOGUILineResizer(groupBoxHisto, "line1");	
	line1->SetMinWidth(200);
	line1->setGeometry( QRect( 263, 8, 7, 131 ) );
	connect(line1, SIGNAL(resizeInfoElem(int)), this, SLOT(changeSizeInfoFirst(int)));

	line2 = new TMOGUILineResizer(groupBoxStat, "line2");	
	line2->SetMinWidth(113);
	line2->setGeometry( QRect( 125, 8, 7, 131 ) );
	connect(line2, SIGNAL(resizeInfoElem(int)), this, SLOT(changeSizeInfoSecond(int)));

	line3 = new TMOGUILineResizer(groupBoxLocalTool, "line3");	
	line3->SetMinWidth(405);
	line3->setGeometry( QRect( 505, 8, 7, 131 ) );
	connect(line3, SIGNAL(resizeInfoElem(int)), this, SLOT(changeSizeInfoThird(int)));

    languageChange();
    resize( QSize(949, 160).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );	
}

/*
 *  Destroys the object and frees any allocated resources
 */
TMOGUIInfoPanel::~TMOGUIInfoPanel()
{
}

void TMOGUIInfoPanel::languageChange()
{
    setCaption( tr( "InfoPanel" ) ); 

	groupBoxHisto->setTitle( tr( "Histogram" ) );
		textLabelBlack->setText( tr( "Black:" ) );
		textLabelWhite->setText( tr( "White:" ) );
		textLabelGammma->setText( tr( "Gamma:" ) );
		textLabelAvg->setText( tr( "Average:" ) );
		pBlack->setText( QString::null );
		pWhite->setText( QString::null );
		pGamma->setText( QString::null );
		pAverage->setText( QString::null );

    groupBoxStat->setTitle( tr( "Statistics" ) );
		textLabelAbove->setText( tr( "Above:" ) );
		textLabelBelow->setText( tr( "Below:" ) );
		textLabelVisible->setText( tr( "Visible:" ) );
		textLabelExtreme->setText( tr( "Extreme:" ) );
		pAbove->setText( QString::null );
		pBelow->setText( QString::null );
		pVisible->setText( QString::null );
		pExtreme->setText( QString::null );

    groupBoxLocalTool->setTitle( tr( "Local info tool statistics" ) );
		textLabelNeighb->setText( tr( "Neighbourhood" ) );		
		textLabelAvgLum->setText( tr( "Avg. Luminance:" ) );
		textLabelMinLum->setText( tr( "Min. Luminace:" ) );
		textLabelMaxLum->setText( tr( "Max. Luminance:" ) );
		textLabelAvgCol->setText( tr( "Avg. Color:" ) );
		pAvgLum->setText( QString::null );
		pMaxLum->setText( QString::null );
		pMinLum->setText( QString::null );
		pAvgCol->setText( QString::null );	
    
		ptextLabelCursor->setText( tr( "Cursor" ) );
		ptextLabelLum->setText( tr( "Luminance:" ) );
		ptextLabelColor->setText( tr( "Color:" ) );
		pLum->setText( QString::null );
		pColor->setText( QString::null );
}

void TMOGUIInfoPanel::changeSizeInfoFirst(int size)
{	
	groupBoxHisto->setGeometry(groupBoxHisto->x(), groupBoxHisto->y(), groupBoxHisto->width() + size, groupBoxHisto->height());
	pBlack->setGeometry(pBlack->x(), pBlack->y(), pBlack->width() + size, pBlack->height());
	pWhite->setGeometry(pWhite->x(), pWhite->y(), pWhite->width() + size, pWhite->height());
	pGamma->setGeometry(pGamma->x(), pGamma->y(), pGamma->width() + size, pGamma->height());
    pAverage->setGeometry(pAverage->x(), pAverage->y(), pAverage->width() + size, pAverage->height());
	groupBoxStat->move(groupBoxStat->x() + size, groupBoxStat->y());
	groupBoxLocalTool->move(groupBoxLocalTool->x() + size, groupBoxLocalTool->y());
	this->setGeometry(this->x(), this->y(), this->width() + size, this->height());
}

void TMOGUIInfoPanel::changeSizeInfoSecond(int size)
{	
	groupBoxStat->setGeometry(groupBoxStat->x(), groupBoxStat->y(), groupBoxStat->width() + size, groupBoxStat->height());
	pBelow->setGeometry(pBelow->x(), pBelow->y(), pBelow->width() + size, pBelow->height());
	pAbove->setGeometry(pAbove->x(), pAbove->y(), pAbove->width() + size, pAbove->height());
	pVisible->setGeometry(pVisible->x(), pVisible->y(), pVisible->width() + size, pVisible->height());
    pExtreme->setGeometry(pExtreme->x(), pExtreme->y(), pExtreme->width() + size, pExtreme->height());
	groupBoxLocalTool->move(groupBoxLocalTool->x() + size, groupBoxLocalTool->y());
	this->setGeometry(this->x(), this->y(), this->width() + size, this->height());
}

void TMOGUIInfoPanel::changeSizeInfoThird(int size)
{	
	int sizeHalf = size / 2;
	groupBoxLocalTool->setGeometry(groupBoxLocalTool->x(), groupBoxLocalTool->y(), groupBoxLocalTool->width() + size, groupBoxLocalTool->height());
	pAvgCol->setGeometry(pAvgCol->x(), pAvgCol->y(), pAvgCol->width() + sizeHalf, pAvgCol->height());
	pMinLum->setGeometry(pMinLum->x(), pMinLum->y(), pMinLum->width() + sizeHalf, pMinLum->height());
	pMaxLum->setGeometry(pMaxLum->x(), pMaxLum->y(), pMaxLum->width() + sizeHalf, pMaxLum->height());
	pAvgLum->setGeometry(pAvgLum->x(), pAvgLum->y(), pAvgLum->width() + sizeHalf, pAvgLum->height());
	pLum->setGeometry(pLum->x(), pLum->y(), pLum->width() + sizeHalf, pLum->height());
	pLum->move(pLum->x() + sizeHalf, pLum->y());
    pColor->setGeometry(pColor->x(), pColor->y(), pColor->width() + sizeHalf, pColor->height());
	pColor->move(pColor->x() + sizeHalf, pColor->y());
	ptextLabelColor->move(ptextLabelColor->x() + sizeHalf, ptextLabelColor->y());
	ptextLabelLum->move(ptextLabelLum->x() + sizeHalf, ptextLabelLum->y());
	ptextLabelCursor->move(ptextLabelCursor->x() + sizeHalf, ptextLabelCursor->y());
	this->setGeometry(this->x(), this->y(), this->width() + size, this->height());
}