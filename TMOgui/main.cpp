#include <QApplication>
#include <QStyleFactory>
#include "TMOGUIStyle.h"
#include "tmoguiwindow.h"

int main( int argc, char** argv )
{

    Q_INIT_RESOURCE (icons);

    QApplication::setStyle(new TMOGUIStyle);

    QApplication app( argc, argv );
    app.setAttribute(Qt::AA_DontShowIconsInMenus, false);

    TMOGUIWindow w;


//  app.setStyle(QStyleFactory::create("Fusion"));
//    app.setStyle(new TMOGUIStyle);
//	app.setMainWidget(&window);
    w.show();

        if(argc>=2)w.openFile(argv[1]);
    return app.exec();
}

