#include <iostream>
#include <QApplication>
#include <QTextEdit>

using namespace std;
int main(int argc, char *argv[])
{
   QApplication app(argc, argv);
	QTextEdit textEdit;
	textEdit.show();
    app.exec();
  
	std::cout << " Finished " << endl;
	return 0;
}
