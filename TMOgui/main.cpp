#include <QApplication>
#include <QStyleFactory>
#include "tmoguiwindow.h"

int main( int argc, char** argv )
{
	QApplication app( argc, argv );

    TMOGUIWindow w;


//  app.setStyle(QStyleFactory::create("Fusion"));
//	app.setMainWidget(&window);

    w.show();

        if(argc>=2)w.openFile(argv[1]);
	return app.exec();
}

