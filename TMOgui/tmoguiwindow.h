#include <qmainwindow.h>
#include <qlist.h>
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

class QWorkspace;
class TMOGUIResource;
class QSplitter;
class TMOGUITransformation;
class TMOGUIImageSize;
class QAssistantClient;


class TMOGUIWindow : public QMainWindow
{
	Q_OBJECT
public:
	TMOGUIWindow( QWidget* parent = 0, const char* name = 0, WFlags f = WType_TopLevel );
	virtual ~TMOGUIWindow();
	void openFile(QString fileName);
	
protected:
	virtual TMOGUIImage* GetActiveImage();
	virtual TMOGUIImage* GetNewImage(const QString &name);
	virtual TMOGUIImage* FindImage(QString name);
	virtual int SavePosition();
	virtual int LoadPosition();
	virtual int Create();
	virtual void WindowChangedToolActivated(TMOGUIImage * pImage);
	QWorkspace *pWorkspace;
	TMOGUIMenu *pMenu;
	TMOGUIZoomTool *pTools;
	TMOGUIInfoToolBar *pInfoTool;
	TMOGUIFileToolBar *pFileTool;
	TMOGUIStatus *pStatus;
	TMOGUIProgressBar *pProgress;
	TMOGUIRightBar *pRight;
	TMOGUIInfo *pInfo;
	TMOGUIImage *pImages[3];
	TMOGUIImageSize *pDialog;
	TMOGUIInfoTool* iTool;

	QList<TMOGUIImage> listImage;
	int iFlags;
	double dRatio;
	double pColors[6];
	double pMargins[4];
	int iWidth, iHeight, iOperation;
	QString sFileName, sPrevFileName;
	bool bZoom;
	QSplitter *pSplitter;
	QSplitter *pRightSplitter;
	QAssistantClient *assistant;

	
	
public slots:
	void newFile();
	void exitFile();
	void openFile();
	void openFile(int);
	void closeFile();
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
	void windowChanged(QWidget*);
	void ImageSizeConstrain(bool bChecked);
	void NewImageConstant(bool bSelected);
	void ImageSizeWidth(const QString &s);
	void ImageSizeHeight(const QString &s);
	void extractLumCommand();
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
		
signals:
	void imageSelected(TMOGUIImage*);
	void signalImageSizeWidth(const QString&);
	void signalImageSizeHeight(const QString&);
	void signalMergeComRed(const QPixmap&);
	void signalMergeComGreen(const QPixmap&);
	void signalMergeComBlue(const QPixmap&);
};

