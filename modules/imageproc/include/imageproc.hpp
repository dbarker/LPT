/**
 * @file imageproc.hpp
 * Image processing module declaration
 */

#ifndef IMAGEPROC_H_
#define IMAGEPROC_H_

#include "core.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
namespace lpt {

using namespace std;

class ImageProcess;
class ImageProcessor;
class SubtractBackground;
class Threshold;
class Erode;
class GaussianBlur;

class Detector;
class FindContoursDetector;
class GoodFeaturesToTrackDetector;

/**
 * @brief Calculates undistorted particle coordinates using OpenCV function undistorPoints
 *
 * @param camera The camera to obtain coefficients
 * @param frame The frame to be processed
 */
void undistortPoints(const lpt::Camera& camera, lpt::ImageFrame& frame);

/**
 * @brief The ImageProcess class
 * Base class for various image processing operations
 */
class ImageProcess {
public:
    typedef std::shared_ptr<ImageProcess> Ptr;

    virtual void addControls() = 0;
    virtual void process(cv::Mat& image) = 0;
    virtual ~ImageProcess();
};

/**
 * @brief The ImageProcessor class
 * The class holds image processes and apply them to
 * the target images
 */
class ImageProcessor {
public:
	typedef std::shared_ptr<ImageProcessor> Ptr;

    /**
     * @brief Creates a shared_ptr to this object
     *
     * @return a shared_ptr to ImageProcessor class
     */
    static inline ImageProcessor::Ptr create() { return std::make_shared<ImageProcessor>(); }

    /**
     * @brief ImageProcessor constructor
     */
    ImageProcessor();

    /**
     * @brief Applies image processes to the target image
     *
     * @param image The image to be processed
     */
    void processImage( cv::Mat& image );

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Adds processes to the image processor
     *
     * @param process The image process to be added
     */
    void addProcess(ImageProcess::Ptr process);

    /**
     * @brief showResults
     *
     * @param window_name
     */
	void showResults( const string& window_name );

private:
    vector<ImageProcess::Ptr> m_processes;
};

/**
 * @brief The SubtractBackground class
 * Derived class from the ImageProcess class
 */
class SubtractBackground : public ImageProcess {
public:
    typedef std::shared_ptr<SubtractBackground> Ptr;
    static inline SubtractBackground::Ptr create() { return make_shared<SubtractBackground>(); }
    /**
     * @brief SubtractBackground construcotr
     */
    SubtractBackground();
    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Subtracts the background
     *
     * @param image The image to be processed
     */
    inline void process(cv::Mat& image) {
        image = image - m_background;
	}

    /**
     * @brief SubtractBackground destructor
     */
    virtual ~SubtractBackground();

private:
    cv::Mat m_background;
};

/**
 * @brief The Threshold class
 * Derived class from the ImageProcess class
 */
class Threshold : public ImageProcess {
public:
    typedef std::shared_ptr<Threshold> Ptr;
    static inline Threshold::Ptr create(int threshold = 100, int max_threshold = 255) {
        return std::make_shared<Threshold>(threshold, max_threshold);
    }

    /**
     * @brief Threshold constructor
     *
     * @param threshold threshold_value
     * @param max max_threshold
     */
    Threshold(int threshold = 100, int max_threshold = 255);

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Applies the threshold
     *
     * @param image The image to be processed
     */
	inline void process(cv::Mat& image) {
        cv::threshold(image, image, m_threshold, m_max_threshold, cv::THRESH_BINARY);
	}

    /**
     * @brief Threshold destroctor
     */
    virtual ~Threshold();

private:
    int m_threshold;
    int m_max_threshold;
};

/**
 * @brief The Erode class
 * Derived class from the ImageProcess class
 */
class Erode : public ImageProcess {
public:
    typedef std::shared_ptr<Erode> Ptr;
    static inline Erode::Ptr create(int iterations = 1, int max_iterations = 5) {
        return std::make_shared<Erode>(iterations, max_iterations);
    }

    /**
     * @brief Erode consturctor
     *
     * @param iterations Number of iterations
     * @param max_iterations Max number of iterations
     */
    Erode(int iterations = 1, int max_iterations = 5);

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Applies Erosion operation
     *
     * @param image The image to be processed
     */
	inline void process(cv::Mat& image) {
        cv::erode(image, image, cv::Mat(), cv::Point(-1,-1), m_iterations);
	}

    /**
     * @brief Erode destructor
     */
    virtual ~Erode();

private:
    int m_iterations;
    int m_max_iterations;
};

/**
 * @brief The EqualizeHistogram class
 * Derived class from the ImageProcess class
 */
class EqualizeHistogram : public ImageProcess {
public:
    typedef std::shared_ptr<EqualizeHistogram> Ptr;
    static inline EqualizeHistogram::Ptr create() {
        return make_shared<EqualizeHistogram>();
    }

    /**
     * @brief EqualizeHistogram constructor
     */
    EqualizeHistogram();

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Applies EqualizeHistogram operation to the image
     *
     * @param image The image to be processed
     */
	inline void process(cv::Mat& image){
		cv::equalizeHist(image, image);
	}

    /**
     * @brief EqualizeHistogram destructor
     */
    virtual ~EqualizeHistogram();
};

/**
 * @brief The Dilate class
 * Derived class from ImageProcess
 */
class Dilate : public ImageProcess {
public:
    typedef std::shared_ptr<Dilate> Ptr;
    static inline Dilate::Ptr create(int iterations = 1, int max_iterations = 5) {
        return make_shared<Dilate>(iterations, max_iterations);
    }

    /**
     * @brief Dilate constructor
     *
     * @param iterations Number of iterations
     * @param max_iterations Max number of iterations
     */
    Dilate(int iterations = 1, int max_iterations = 5);

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Applies dilation operation to the image
     *
     * @param image The image to be processed
     */
	inline void process(cv::Mat& image) {
        cv::dilate(image, image, cv::Mat(), cv::Point(-1,-1), m_iterations);
	}

    /**
     * @brief Dilate destructor
     */
    virtual ~Dilate();

private:
    int m_iterations;
    int m_max_iterations;
};

/**
 * @brief The GaussianBlur class
 * Derived class from the ImageProcess class
 */
class GaussianBlur : public ImageProcess {
public:
    typedef std::shared_ptr<GaussianBlur> Ptr;
    static inline GaussianBlur::Ptr create(int kernel_size = 5, double sigma1 = 0.0, double sigma2 = 0.0, int boarder = 4) {
        return std::make_shared<GaussianBlur>(kernel_size, sigma1, sigma2, boarder);
    }

    /**
     * @brief GaussianBlur constructor
     *
     * @param kernal_size
     * @param sigma1
     * @param sigma2
     * @param boarder
     */
    GaussianBlur(int kernel_size = 5, double sigma1 = 0.0, double sigma2 = 0.0, int boarder = 4);

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Applies Gaussian blur to the image
     *
     * @param image The image to be processed
     */
	inline void process(cv::Mat& image) {
        cv::GaussianBlur( image, image, cv::Size(m_kernel_size, m_kernel_size), m_sigma1, m_sigma2, m_boarder_type);
	}

    /**
     * @brief GaussianBlur destructor
     */
    virtual ~GaussianBlur();

private:
    int m_kernel_size;
    double m_sigma1;
    double m_sigma2;
    int m_boarder_type;
};

/**
 * @brief The Detector class
 * Base class for determining particle centroid
 */
class Detector {
public:
	typedef std::shared_ptr<Detector> Ptr;

    virtual void detectFeatures( const cv::Mat& image, vector<lpt::ParticleImage::Ptr>& features, vector<vector<cv::Point>>& contours ) = 0;
	virtual void addControls() = 0;

    /**
     * @brief Draws circles around detected particles
     *
     * @param frame
     */
    void drawResult(ImageFrame& frame);

	void drawContours( cv::Mat& image, vector<vector<cv::Point>> contours );

    /**
     * @brief Detector destructor
     */
    virtual ~Detector();
};

/**
 * @brief The FindContoursDetector class
 * Derived class from the Detector class
 * An intensity-based weighted average method
 */
class FindContoursDetector : public Detector {
public:
	typedef std::shared_ptr<FindContoursDetector> Ptr;

    /**
     * @brief Creates a shared_ptr to this object
     *
     * @return A shared_ptr to the FindContoursDetector class
     */
	static inline FindContoursDetector::Ptr create() { return std::make_shared<lpt::FindContoursDetector>(); }

    /**
     * @brief The Parameters class
     * Specifies the parameters for centroid detection
     */
	class Parameters {
	public:
		Parameters(int mode = cv::RETR_EXTERNAL, int method = cv::CHAIN_APPROX_NONE, int min_area = 2, int max_area = 200 )
			: mode(mode), method(method), min_contour_area(min_area), max_contour_area(max_area) {}
		int mode;
		int method;
		int min_contour_area;
		int max_contour_area;
	} params;
	
    /**
     * @brief FindContoursDetector constructor
     */
    FindContoursDetector();

    /**
     * @brief Detects features in the image
     * Finds contours of particle image using OpenCV function findContours
     * Then calculates particle centroid by weighted-averaging
     *
     * @param image The image to be processed
     * @param features The vector holding detected features
     */
    void detectFeatures(const cv::Mat &image, vector<lpt::ParticleImage::Ptr>& features, vector<vector<cv::Point>>& contours );

	//void detectFeatures(const cv::Mat &image, vector<lpt::ParticleImage::Ptr>& features, vector<vector<cv::Point>>& contours);

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief Draws detected contours
     *
     * @param result_image The image with contours drawn
     * @param contours The contours to be drawn
     */
    //void drawContours(cv::Mat& result_image, vector<vector<cv::Point> > contours);

    /**
     * @brief FindContoursDetector destructor
     */
    virtual ~FindContoursDetector();
};

/**
 * @brief The GoodFeaturesToTrackDetector class
 * Derived class from the Detector class
 * Good feature method from Shindler 2010
 */
class GoodFeaturesToTrackDetector : public Detector {
public:
	typedef std::shared_ptr<GoodFeaturesToTrackDetector> Ptr;

    /**
     * @brief Creates a shared_ptr to this object
     * @return A shared_ptr to the GoodFeaturesToTrackDetector class
     */
	static inline GoodFeaturesToTrackDetector::Ptr create() { return std::make_shared<lpt::GoodFeaturesToTrackDetector>(); }

    /**
     * @brief The Parameters class
     * Specifies the parameters for the Good Feature Detection method
     */
	class Parameters {
    public:
        Parameters ( int max_corners = 10000, double quality = 0.01,
                double min_distance = 5.0, cv::Mat mask = cv::Mat(), int blocksize = 6,
                bool use_harris = false, double k = 0.04 )
				: max_corners(max_corners), quality_level(quality), min_distance(min_distance),
				mask(mask), neighborhood_size(blocksize), use_harris(use_harris), k(k) {}
        int max_corners;
        double quality_level;
        double min_distance;
        cv::Mat mask;
		int neighborhood_size;
        bool use_harris;
        double k;
    } params;

    /**
     * @brief GoodFeaturesToTrackDetector constructor
     */
    GoodFeaturesToTrackDetector();
	
    /**
     * @brief Detects features in the image
     * Finds contours of particle image using OpenCV function findContours
     * Then calculates particle centroid by weighted-averagin
     *
     * @param image The image to be processed
     * @param features The vector holding detected features
     */
    void detectFeatures( const cv::Mat& image, vector<lpt::ParticleImage::Ptr>& features, vector<vector<cv::Point>>& contours );

    /**
     * @brief addControls
     */
    void addControls();

    /**
     * @brief GoodFeaturesToTrackDetector destructor
     */
    virtual ~GoodFeaturesToTrackDetector();
};

/**
 * @brief processImages
 * @param camera
 * @param processor
 * @param detector
 */
void processImages(
		Camera& camera,
		ImageProcessor& processor,
		Detector& detector);

/**
 * @brief testDetectedParticles
 * @param true_particles
 * @param detected_particles
 */
void testDetectedParticles(
		vector<ParticleImage::Ptr>& true_particles,
		vector<ParticleImage::Ptr>& detected_particles);

}/* NAMESPACE_PT */

#endif /*IMAGEPROC_H_*/

