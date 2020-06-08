#include <QMainWindow>
#include <qlist.h>
//Added by qt3to4:
#include <QPixmap>
#include "TMOGUIMenu.h"
#include "TMOGUIZoomTool.h"
#include "TMOGUIStatus.h"
#include "TMOGUIImage.h"
#include "TMOGUIRightBar.h"
#include "TMOGUIFilters.h"
#include "TMOGUIInfo.h"
#include "TMOGUIInfoToolBar.h"
#include "TMOGUIFileToolBar.h"
#include "TMOGUIProgressBar.h"
#include "ui_resources.h"
#include "ui_resources1.h"
#include "ui_resources2.h"
#include "ui_resources3.h"
#include "ui_resources4.h"
#include "ui_resources5.h"
#include "ui_resources6.h"

class QMdiArea;
class QMdiSubWindow;
class TMOGUIResource;
class QSplitter;
class TMOGUITransformation;
// class QAssistantClient;


class TMOGUIWindow : public QMainWindow
{
	Q_OBJECT
public:
    TMOGUIWindow( QWidget* parent = 0, const char* name = 0, Qt::WindowFlags f = 0 );
	virtual ~TMOGUIWindow();
    //void openFile(QString fileName);


protected:
	virtual TMOGUIImage* GetActiveImage();
	virtual TMOGUIImage* GetNewImage(const QString &name);
	virtual TMOGUIImage* FindImage(QString name);
	virtual int SavePosition();
	virtual int LoadPosition();
	virtual int Create();
	virtual void WindowChangedToolActivated(TMOGUIImage * pImage);
    void mySetItemChecked(TMOGUIMenu *menu, QMenu *item, int id, bool check);
    QMdiArea *pWorkspace;
    TMOGUIMenu *pMenu;
    TMOGUIZoomTool *pTools;
    TMOGUIInfoToolBar *pInfoTool;
    TMOGUIFileToolBar *pFileTool;
    TMOGUIStatus *pStatus;
    TMOGUIProgressBar *pProgress;
    TMOGUIRightBar *pRight;
    TMOGUIInfo *pInfo;
    TMOGUIImage *pImages[3];
    Ui::TMOGUIImageSize *pDialog;
    TMOGUIInfoTool* iTool;

    QList<TMOGUIImage*> listImage;
	int iFlags;
	double dRatio;
	double pColors[6];
	double pMargins[4];
	int iWidth, iHeight, iOperation;
	QString sFileName, sPrevFileName;
	bool bZoom;
	QSplitter *pSplitter;
	QSplitter *pRightSplitter;
    //QAssistantClient *assistant;

private:
    void setup();
	
	
public slots:
	void newFile();
	void exitFile();
	void openFile();
	void openFile(int);
    void openFile(QString);
	void closeFile();
    void closeActiveWindow();
	void saveFile();
	void saveasFile();
	void saveallFile();
	void closeallWindow();
	void activateWindow(int id);
	void undoEdit();
	void transform();
	void viewHistogram();
	void viewRight();
	void viewInfo();
	void duplicateCommand();
	void sizeCommand();
    void windowChanged(QMdiSubWindow* pWidget);
	void ImageSizeConstrain(bool bChecked);
	void NewImageConstant(bool bSelected);
	void ImageSizeWidth(const QString &s);
	void ImageSizeHeight(const QString &s);
	void extractLumCommand();
    void extractRed();
    void extractGreen();
    void extractBlue();
	void extractComCommand(int);
	void mergeCommand();
	void operationCommand();
	void OperationFirst(int iImage);
	void OperationSecond(int iImage);
	void MergeComponentsRed(int iImage);
	void MergeComponentsGreen(int iImage);
	void MergeComponentsBlue(int iImage);
	void ImageOperation(int iOperation);
	void SetRMin(const QString&);
	void SetGMin(const QString&);
	void SetBMin(const QString&);
	void SetRMax(const QString&);
	void SetGMax(const QString&);
	void SetBMax(const QString&);
	void SetFileName(const QString&);
	void SetOperation(int iOp);
	void SetWidth(const QString& s);
	void SetHeight(const QString& s);
	void pageFile();
	void printFile();
	void finishTransform();
	void zoom();
	void fitToScreen();
	void zoomIn();
	void zoomOut();
	void fitToWidth();
	void fitToHeight();
	void activateInfoTool(bool on);
	void showToolSetting();
	void setAssistantArguments();
    void showAssistantErrors( const QString &err );    
    void showHelp();
    void refreshWindowsList();
		
signals:
	void imageSelected(TMOGUIImage*);
	void signalImageSizeWidth(const QString&);
	void signalImageSizeHeight(const QString&);
	void signalMergeComRed(const QPixmap&);
	void signalMergeComGreen(const QPixmap&);
	void signalMergeComBlue(const QPixmap&);
};

