package bookexamples;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.util.List;

import opennlp.tools.cmdline.parser.ParserTool;
import opennlp.tools.parser.Parse;
import opennlp.tools.parser.Parser;
import opennlp.tools.parser.ParserFactory;
import opennlp.tools.parser.ParserModel;
import edu.stanford.nlp.ling.CoreLabel;
import edu.stanford.nlp.ling.Sentence;
import edu.stanford.nlp.parser.lexparser.LexicalizedParser;
import edu.stanford.nlp.process.CoreLabelTokenFactory;
import edu.stanford.nlp.process.PTBTokenizer;
import edu.stanford.nlp.process.Tokenizer;
import edu.stanford.nlp.process.TokenizerFactory;
import edu.stanford.nlp.trees.GrammaticalStructure;
import edu.stanford.nlp.trees.GrammaticalStructureFactory;
import edu.stanford.nlp.trees.Tree;
import edu.stanford.nlp.trees.TreebankLanguagePack;
import edu.stanford.nlp.trees.TypedDependency;

/**
 * Chapter 7 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter7 {

	public static void opennlpParser() {
		String fileLocation = "add file directory/en-parser-chunking.bin";
		try (InputStream modelInputStream = new FileInputStream(fileLocation);) {
			ParserModel model = new ParserModel(modelInputStream);
			Parser parser = ParserFactory.create(model);
			String sentence = "The cow jumped over the moon";
			Parse parses[] = ParserTool.parseLine(sentence, parser, 3);
			for (Parse parse : parses) {
				parse.show();
				Parse children[] = parse.getChildren();
				for (Parse parseElement : children) {
					System.out.println(parseElement.getText());
					System.out.println(parseElement.getType());
					Parse tags[] = parseElement.getTagNodes();
					System.out.println("Tags");
					for (Parse tag : tags) {
						System.out.println("[" + tag + "] type: "
								+ tag.getType() + " Probability "
								+ tag.getProb() + " Label: " + tag.getLabel());
					}
				}
			}
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void stanfordParser() {
		String parserModel = "add file directory/englishPCFG.ser.gz";
		LexicalizedParser lexicalizedParser = LexicalizedParser
				.loadModel(parserModel);
		String sentenceArray[] = { "The", "cow", "jumped", "over", "the",
				"moon", "." };
		List<CoreLabel> words = Sentence.toCoreLabelList(sentenceArray);
		Tree parseTree = lexicalizedParser.apply(words);
		parseTree.pennPrint();
		System.out.println();
		String sentence = "The cow jumped over the moon.";
		TokenizerFactory<CoreLabel> tokenizerFactory = PTBTokenizer.factory(
				new CoreLabelTokenFactory(), "");
		Tokenizer<CoreLabel> tokenizer = tokenizerFactory
				.getTokenizer(new StringReader(sentence));
		List<CoreLabel> wordList = tokenizer.tokenize();
		parseTree = lexicalizedParser.apply(wordList);
		TreebankLanguagePack tlp = lexicalizedParser.treebankLanguagePack();
		GrammaticalStructureFactory gsf = tlp.grammaticalStructureFactory();
		GrammaticalStructure gs = gsf.newGrammaticalStructure(parseTree);
		List<TypedDependency> tdl = gs.typedDependenciesCCprocessed();
		System.out.println(tdl);
		for (TypedDependency dependency : tdl) {
			System.out.println("Governor Word: [" + dependency.gov()
					+ "] Relation: [" + dependency.reln().getLongName()
					+ "] Dependent Word: [" + dependency.dep() + "]");
		}
		System.out.println();
	}

	public static void stanfordQuestionAnswer() {
		String question = "Who is the 32nd president of the United States?";
		String parserModel = "add file directory/englishPCFG.ser.gz";
		LexicalizedParser lexicalizedParser = LexicalizedParser
				.loadModel(parserModel);
		TokenizerFactory<CoreLabel> tokenizerFactory = PTBTokenizer.factory(
				new CoreLabelTokenFactory(), "");
		Tokenizer<CoreLabel> tokenizer = tokenizerFactory
				.getTokenizer(new StringReader(question));
		List<CoreLabel> wordList = tokenizer.tokenize();
		Tree parseTree = lexicalizedParser.apply(wordList);
		TreebankLanguagePack tlp = lexicalizedParser.treebankLanguagePack();
		GrammaticalStructureFactory gsf = tlp.grammaticalStructureFactory();
		GrammaticalStructure gs = gsf.newGrammaticalStructure(parseTree);
		List<TypedDependency> tdl = gs.typedDependenciesCCprocessed();
		System.out.println(tdl);
		for (TypedDependency dependency : tdl) {
			System.out.println("Governor Word: [" + dependency.gov()
					+ "] Relation: [" + dependency.reln().getLongName()
					+ "] Dependent Word: [" + dependency.dep() + "]");
		}
		System.out.println();
	}

	/**
	 * Runs chapter 7 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		opennlpParser();
		stanfordParser();
		stanfordQuestionAnswer();
	}

}