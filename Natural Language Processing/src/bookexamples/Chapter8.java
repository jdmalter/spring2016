package bookexamples;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Properties;
import java.util.Set;

import opennlp.tools.sentdetect.SentenceDetectorME;
import opennlp.tools.sentdetect.SentenceModel;
import opennlp.tools.tokenize.WhitespaceTokenizer;

import org.apache.xmlbeans.XmlException;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import de.l3s.boilerpipe.BoilerpipeProcessingException;
import de.l3s.boilerpipe.document.TextDocument;
import de.l3s.boilerpipe.sax.BoilerpipeSAXInput;
import de.l3s.boilerpipe.sax.HTMLDocument;
import de.l3s.boilerpipe.sax.HTMLFetcher;
import edu.stanford.nlp.pipeline.Annotation;
import edu.stanford.nlp.pipeline.StanfordCoreNLP;

/**
 * Chapter 8 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter8 {

	public static void boilerpipeHTML() {
		try {
			URL url = new URL("https://en.wikipedia.org/wiki/Berlin");
			HTMLDocument htmlDoc = HTMLFetcher.fetch(url);
			InputSource is = htmlDoc.toInputSource();
			TextDocument document = new BoilerpipeSAXInput(is)
					.getTextDocument();
			System.out.println(document.getText(true, true));
		} catch (MalformedURLException ex) {
			System.out.println("MalformedURL");
		} catch (BoilerpipeProcessingException | SAXException | IOException ex) {
			System.out.println("Other");
		} catch (NoClassDefFoundError er) {
			System.out.println("NoClassDefFound");
		}
	}

	public static void apachePOI() {
		try {
			FileInputStream fis = new

			FileInputStream("TestDocument.docx");
			POITextExtractor textExtractor = ExtractorFactory
					.createExtractor(fis);
			System.out.println(textExtractor.getText());
			POITextExtractor metaExtractor = textExtractor
					.getMetadataTextExtractor();
			System.out.println(metaExtractor.getText());
			CoreProperties coreProperties = properties.getCoreProperties();
			System.out.println(properties.getCorePropertiesText());
			ExtendedProperties extendedProperties = properties
					.getExtendedProperties();
			System.out.println(properties.getExtendedPropertiesText());
		} catch (IOException ex) {
			System.out.println("IO");
		} catch (OpenXML4JException | XmlException ex) {
			System.out.println("Other");
		}
	}

	public static void apachePDF() {
		try {
			File file = new File("TestDocument.pdf");
			PDDocument pdDocument = PDDocument.load(file);
			PDFTextStripper stripper = new PDFTextStripper();
			String text = stripper.getText(pdDocument);
			System.out.println(text);
			pdDocument.close();
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void stanfordPipeline() {
		String text = "The robber took the cash and ran.";
		Properties props = new Properties();
		props.put("annotators",
				"tokenize, ssplit, pos, lemma, ner, parse, dcoref");
		StanfordCoreNLP pipeline = new StanfordCoreNLP(props);
	}

	public static void opennlpPipeline() {
		try (InputStream is = new FileInputStream(new File(
				"add file directory/en-sent.bin"));
				FileReader fr = new FileReader("Twenty Thousands.txt");
				BufferedReader br = new BufferedReader(fr)) {
			SentenceModel model = new SentenceModel(is);
			SentenceDetectorME detector = new SentenceDetectorME(model);
			String line;
			StringBuilder sb = new StringBuilder();
			while ((line = br.readLine()) != null) {
				sb.append(line + " ");
			}
			String sentences[] = detector.sentDetect(sb.toString());
			for (int i = 0; i < sentences.length; i++) {
				sentences[i] = sentences[i].toLowerCase();
			}

		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	/**
	 * Runs chapter 8 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		boilerpipeHTML();
		apachePOI();
		apachePDF();
		stanfordPipeline();
		opennlpPipeline();
	}

}