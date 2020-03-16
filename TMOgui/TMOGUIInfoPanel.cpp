
#include "TMOGUIInfoPanel.h"
#include "TMOGUILineResizer.h"

#include <qlabel.h>
#include <QGroupBox>
#include <QGridLayout>
#include <qtooltip.h>


TMOGUIInfoPanel::TMOGUIInfoPanel( QWidget* parent, const char* name )
    : QWidget( parent)
{
    groupBoxStat = new QGroupBox( this );//, "groupBoxStat" );
    groupBoxStat->setGeometry( QRect( 290, 10, 130, 140 ) );
    QGridLayout *gridStat = new QGridLayout(groupBoxStat);


    pAbove = new QLabel( groupBoxStat);//, "pAbove" );
    pAbove->setGeometry( QRect( 70, 50, 49, 20 ) );


    textLabelBelow = new QLabel( groupBoxStat);//, "textLabelBelow" );
    textLabelBelow->setGeometry( QRect( 20, 30, 49, 20 ) );
    QFont textLabelBelow_font(  textLabelBelow->font() );
    textLabelBelow_font.setBold( true );
    textLabelBelow->setFont( textLabelBelow_font ); 


    textLabelAbove = new QLabel( groupBoxStat);//, "textLabelAbove" );
    textLabelAbove->setGeometry( QRect( 20, 50, 49, 20 ) );
    QFont textLabelAbove_font(  textLabelAbove->font() );
    textLabelAbove_font.setBold( true );
    textLabelAbove->setFont( textLabelAbove_font );


    textLabelVisible = new QLabel( groupBoxStat);//, "textLabelVisible" );
    textLabelVisible->setGeometry( QRect( 20, 70, 49, 20 ) );
    QFont textLabelVisible_font(  textLabelVisible->font() );
    textLabelVisible_font.setBold( true );
    textLabelVisible->setFont( textLabelVisible_font ); 


    textLabelExtreme = new QLabel( groupBoxStat);//, "textLabelExtreme" );
    textLabelExtreme->setGeometry( QRect( 20, 100, 49, 20 ) );
    QFont textLabelExtreme_font(  textLabelExtreme->font() );
    textLabelExtreme_font.setBold( true );
    textLabelExtreme->setFont( textLabelExtreme_font ); 

    pExtreme = new QLabel( groupBoxStat);//, "pExtreme" );
    pExtreme->setGeometry( QRect( 70, 100, 49, 20 ) );

    pVisible = new QLabel( groupBoxStat);//, "pVisible" );
    pVisible->setGeometry( QRect( 70, 70, 49, 20 ) );

    pBelow = new QLabel( groupBoxStat);//, "pBelow" );
    pBelow->setGeometry( QRect( 70, 30, 49, 20 ) );

    gridStat->addWidget(textLabelBelow, 0, 0);
    gridStat->addWidget(textLabelAbove, 1, 0);
    gridStat->addWidget(textLabelVisible, 2, 0);
    gridStat->addWidget(textLabelExtreme, 3, 0);

    gridStat->addWidget(pBelow, 0, 1);
    gridStat->addWidget(pAbove, 1, 1);
    gridStat->addWidget(pVisible, 2, 1);
    gridStat->addWidget(pExtreme, 3, 1);

    //----------------------------------------------------------------

    groupBoxHisto = new QGroupBox( this);//, "groupBoxHisto" );
    groupBoxHisto->setGeometry( QRect( 10, 10, 270, 140 ) );
    QGridLayout *gridHisto = new QGridLayout(groupBoxStat);

    pGamma = new QLabel( groupBoxHisto );//, "pGamma" );
    pGamma->setGeometry( QRect( 80, 70, 180, 20 ) );

    pAverage = new QLabel( groupBoxHisto );//, "pAverage" );
    pAverage->setGeometry( QRect( 80, 100, 180, 20 ) );

    pWhite = new QLabel( groupBoxHisto );//, "pWhite" );
    pWhite->setGeometry( QRect( 80, 50, 180, 20 ) );

    pBlack = new QLabel( groupBoxHisto );//, "pBlack" );
    pBlack->setGeometry( QRect( 80, 30, 180, 20 ) );

    textLabelBlack = new QLabel( groupBoxHisto );//, "textLabelBlack" );
    textLabelBlack->setGeometry( QRect( 20, 30, 49, 20 ) );
    QFont textLabelBlack_font(  textLabelBlack->font() );
    textLabelBlack_font.setBold( true );
    textLabelBlack->setFont( textLabelBlack_font ); 

    textLabelWhite = new QLabel( groupBoxHisto );//, "textLabelWhite" );
    textLabelWhite->setGeometry( QRect( 20, 50, 49, 20 ) );
    QFont textLabelWhite_font(  textLabelWhite->font() );
    textLabelWhite_font.setBold( true );
    textLabelWhite->setFont( textLabelWhite_font ); 

    textLabelGammma = new QLabel( groupBoxHisto );//, "textLabelGammma" );
    textLabelGammma->setGeometry( QRect( 20, 70, 49, 20 ) );
    QFont textLabelGammma_font(  textLabelGammma->font() );
    textLabelGammma_font.setBold( true );
    textLabelGammma->setFont( textLabelGammma_font ); 

    textLabelAvg = new QLabel( groupBoxHisto );//, "textLabelAvg" );
    textLabelAvg->setGeometry( QRect( 20, 100, 49, 20 ) );
    QFont textLabelAvg_font(  textLabelAvg->font() );
    textLabelAvg_font.setBold( true );
    textLabelAvg->setFont( textLabelAvg_font ); 

    gridHisto->addWidget(textLabelBlack, 0, 0);
    gridHisto->addWidget(textLabelWhite, 1, 0);
    gridHisto->addWidget(textLabelGammma, 2, 0);
    gridHisto->addWidget(textLabelAvg, 3, 0);

    gridHisto->addWidget(pBlack, 0, 1);
    gridHisto->addWidget(pWhite, 1, 1);
    gridHisto->addWidget(pGamma, 2, 1);
    gridHisto->addWidget(pAverage, 3, 1);

    //----------------------------------------------------------------

    groupBoxLocalTool = new QGroupBox( this);//, "groupBoxLocalTool" );
    groupBoxLocalTool->setGeometry( QRect( 430, 10, 510, 140 ) );
    QGridLayout *gridLocalTool = new QGridLayout(groupBoxStat);

    textLabelAvgLum = new QLabel( groupBoxHisto );//, "textLabelAvgLum" );
    textLabelAvgLum->setGeometry( QRect( 20, 50, 95, 20 ) );
    QFont textLabelAvgLum_font(  textLabelAvgLum->font() );
    textLabelAvgLum_font.setBold( true );
    textLabelAvgLum->setFont( textLabelAvgLum_font ); 
    textLabelAvgLum->setToolTip("Average value of luminance in the area");

    pMaxLum = new QLabel( groupBoxHisto );//, "pMaxLum" );
    pMaxLum->setGeometry( QRect( 120, 70, 150, 20 ) );

    pMinLum = new QLabel( groupBoxHisto );//, "pMinLum" );
    pMinLum->setGeometry( QRect( 120, 90, 150, 20 ) );

    textLabelAvgCol = new QLabel( groupBoxHisto );//, "textLabelAvgCol" );
    textLabelAvgCol->setGeometry( QRect( 20, 110, 70, 20 ) );
    QFont textLabelAvgCol_font(  textLabelAvgCol->font() );
    textLabelAvgCol_font.setBold( true );
    textLabelAvgCol->setFont( textLabelAvgCol_font ); 
    textLabelAvgCol->setToolTip("Average value of color in the area <red, green, blue>");

    textLabelMinLum = new QLabel( groupBoxHisto );//, "textLabelMinLum" );
    textLabelMinLum->setGeometry( QRect( 20, 90, 86, 20 ) );
    QFont textLabelMinLum_font(  textLabelMinLum->font() );
    textLabelMinLum_font.setBold( true );
    textLabelMinLum->setFont( textLabelMinLum_font ); 
    textLabelMinLum->setToolTip("Minimum value of luminance in the area");

    textLabelMaxLum = new QLabel( groupBoxHisto );//, "textLabelMaxLum" );
    textLabelMaxLum->setGeometry( QRect( 20, 70, 96, 20 ) );
    QFont textLabelMaxLum_font(  textLabelMaxLum->font() );
    textLabelMaxLum_font.setBold( true );
    textLabelMaxLum->setFont( textLabelMaxLum_font ); 
    textLabelMaxLum->setToolTip("Maximum value of luminance in the area");

    textLabelNeighb = new QLabel( groupBoxHisto );//, "textLabelNeighb" );
    textLabelNeighb->setGeometry( QRect( 20, 20, 89, 20 ) );
    QFont textLabelNeighb_font(  textLabelNeighb->font() );
    textLabelNeighb_font.setBold( true );
    textLabelNeighb_font.setItalic( true );
    textLabelNeighb_font.setUnderline( true );
    textLabelNeighb->setFont( textLabelNeighb_font );
    textLabelNeighb->setToolTip("Data from the place determined by the shape (circle/square)");

    ptextLabelCursor = new QLabel( groupBoxHisto );//, "ptextLabelCursor" );
    ptextLabelCursor->setGeometry( QRect( 280, 20, 89, 20 ) );
    QFont ptextLabelCursor_font(  ptextLabelCursor->font() );
    ptextLabelCursor_font.setBold( true );
    ptextLabelCursor_font.setItalic( true );
    ptextLabelCursor_font.setUnderline( true );
    ptextLabelCursor->setFont( ptextLabelCursor_font ); 
    ptextLabelCursor->setToolTip("Data from the place below the cursor");

    ptextLabelColor = new QLabel( groupBoxHisto );//, "ptextLabelColor" );
    ptextLabelColor->setGeometry( QRect( 280, 70, 40, 20 ) );
    QFont ptextLabelColor_font(  ptextLabelColor->font() );
    ptextLabelColor_font.setBold( true );
    ptextLabelColor->setFont( ptextLabelColor_font ); 
    ptextLabelColor->setToolTip("Color under the cursor");

    pAvgLum = new QLabel( groupBoxHisto );//, "pAvgLum" );
    pAvgLum->setGeometry( QRect( 120, 50, 150, 20 ) );

    ptextLabelLum = new QLabel( groupBoxHisto );//, "ptextLabelLum" );
    ptextLabelLum->setGeometry( QRect( 280, 50, 70, 20 ) );
    QFont ptextLabelLum_font(  ptextLabelLum->font() );
    ptextLabelLum_font.setBold( true );
    ptextLabelLum->setFont( ptextLabelLum_font ); 
    ptextLabelLum->setToolTip("Luminance under the cursor");

    pLum = new QLabel( groupBoxLocalTool);//, "pLum" );
    pLum->setGeometry( QRect( 350, 50, 150, 20 ) );

    pColor = new QLabel( groupBoxLocalTool);//, "pColor" );
    pColor->setGeometry( QRect( 350, 70, 150, 20 ) );

    pAvgCol = new QLabel( groupBoxLocalTool);//, "pAvgCol" );
    pAvgCol->setGeometry( QRect( 120, 110, 150, 20 ) );

    gridLocalTool->addWidget(textLabelNeighb, 0, 0, 1, 0);
    gridLocalTool->addWidget(textLabelAvgLum, 1, 0);
    gridLocalTool->addWidget(textLabelMinLum, 2, 0);
    gridLocalTool->addWidget(textLabelMaxLum, 3, 0);
    gridLocalTool->addWidget(textLabelAvgCol, 4, 0);

    gridLocalTool->addWidget(ptextLabelCursor, 0, 2, 1, 0);
    gridLocalTool->addWidget(ptextLabelLum, 1, 2);
    gridLocalTool->addWidget(ptextLabelColor, 2, 2);

    //gridLocalTool->addWidget(nullptr, 0, 1);
    gridLocalTool->addWidget(pAvgLum, 1, 1);
    gridLocalTool->addWidget(pMinLum, 2, 1);
    gridLocalTool->addWidget(pMaxLum, 3, 1);
    gridLocalTool->addWidget(pAvgCol, 4, 1);

    gridLocalTool->addWidget(pLum, 1, 3);
    gridLocalTool->addWidget(pColor, 2, 3);


	line1 = new TMOGUILineResizer(groupBoxHisto, "line1");	
    line1->SetMinWidth(250);
	line1->setGeometry( QRect( 263, 8, 7, 131 ) );
    connect(line1, &TMOGUILineResizer::resizeInfoElem, this, &TMOGUIInfoPanel::changeSizeInfoFirst);

	line2 = new TMOGUILineResizer(groupBoxStat, "line2");	
    line2->SetMinWidth(250);
	line2->setGeometry( QRect( 125, 8, 7, 131 ) );
    connect(line2, &TMOGUILineResizer::resizeInfoElem, this, &TMOGUIInfoPanel::changeSizeInfoSecond);

	line3 = new TMOGUILineResizer(groupBoxLocalTool, "line3");	
	line3->SetMinWidth(405);
	line3->setGeometry( QRect( 505, 8, 7, 131 ) );
    connect(line3, &TMOGUILineResizer::resizeInfoElem, this, &TMOGUIInfoPanel::changeSizeInfoThird);


    groupBoxStat->setLayout(gridStat);
    groupBoxHisto->setLayout(gridHisto);
    groupBoxLocalTool->setLayout(gridLocalTool);

    languageChange();
    resize( QSize(949, 160).expandedTo(minimumSizeHint()) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
TMOGUIInfoPanel::~TMOGUIInfoPanel()
{
}

void TMOGUIInfoPanel::languageChange()
{
    setWindowTitle( tr( "InfoPanel" ) );

	groupBoxHisto->setTitle( tr( "Histogram" ) );
		textLabelBlack->setText( tr( "Black:" ) );
		textLabelWhite->setText( tr( "White:" ) );
		textLabelGammma->setText( tr( "Gamma:" ) );
		textLabelAvg->setText( tr( "Average:" ) );
        pBlack->setText( QString() );
        pWhite->setText( QString() );
        pGamma->setText( QString() );
        pAverage->setText( QString() );

    groupBoxStat->setTitle( tr( "Statistics" ) );
		textLabelAbove->setText( tr( "Above:" ) );
		textLabelBelow->setText( tr( "Below:" ) );
		textLabelVisible->setText( tr( "Visible:" ) );
		textLabelExtreme->setText( tr( "Extreme:" ) );
        pAbove->setText( QString() );
        pBelow->setText( QString() );
        pVisible->setText( QString() );
        pExtreme->setText( QString() );

    groupBoxLocalTool->setTitle( tr( "Local info tool statistics" ) );
		textLabelNeighb->setText( tr( "Neighbourhood" ) );		
		textLabelAvgLum->setText( tr( "Avg. Luminance:" ) );
		textLabelMinLum->setText( tr( "Min. Luminace:" ) );
		textLabelMaxLum->setText( tr( "Max. Luminance:" ) );
		textLabelAvgCol->setText( tr( "Avg. Color:" ) );
        pAvgLum->setText( QString() );
        pMaxLum->setText( QString() );
        pMinLum->setText( QString() );
        pAvgCol->setText( QString() );
    
		ptextLabelCursor->setText( tr( "Cursor" ) );
		ptextLabelLum->setText( tr( "Luminance:" ) );
		ptextLabelColor->setText( tr( "Color:" ) );
        pLum->setText( QString() );
        pColor->setText( QString() );
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
