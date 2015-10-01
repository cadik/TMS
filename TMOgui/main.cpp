#include <qapplication.h>
//#include <qmotifplusstyle.h>
#include "tmoguiwindow.h"

int main( int argc, char** argv )
{
	QApplication app( argc, argv );

	TMOGUIWindow window;
//	QApplication::setStyle( new QMotifPlusStyle );
	app.setMainWidget(&window);

	window.show();
        if(argc>=2)window.openFile(argv[1]);
	return app.exec();
}

