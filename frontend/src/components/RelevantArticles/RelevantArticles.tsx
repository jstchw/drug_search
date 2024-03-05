import { useRelevantArticles } from '../../hooks/useRelevantArticles';
import { Books, LinkSimple, X } from '@phosphor-icons/react';
import { Accordion, Button, ButtonGroup, Row, Col, Badge } from 'react-bootstrap';
import { highlightWords } from '../../utils/utils';
import { useUrlParams } from '../../hooks/useUrlParams';
import useArticleStore from '../../stores/articleStore';
import { useEffect } from 'react';
import { capitalizeFirstLetter } from '../../utils/utils';

const RelevantArticles = () => {
  const { relevantArticles, relevantArticlesError } = useRelevantArticles();
  const {
    params: { terms, searchBy },
  } = useUrlParams();

  const articleTerms = useArticleStore((state) => state.articleTerms);

  useEffect(() => {
    useArticleStore.getState().resetArticleTerms();

    terms.forEach((term) => {
      useArticleStore.getState().addArticleTerm(capitalizeFirstLetter(term), searchBy);
    });
  }, [searchBy, terms]);

  if (relevantArticlesError || relevantArticles.length === 0) {
    return;
  }

  return (
    <div className={'mb-5'}>
      <h3 className={'d-flex justify-content-center align-items-center mb-2'}>
        <Books className={'text-secondary'} weight={'light'} />
        <div className={'vr mx-2'} />
        Relevant Articles
      </h3>
      <Row className={'mb-4'}>
        <Col>
          <div className={'text-secondary d-flex justify-content-center mb-2'}>Matching for:</div>
          <div className={'d-flex justify-content-center'}>
            <h5>
              {articleTerms.map((item) => {
                return (
                  <ButtonGroup key={item.term + item.type} className={'mx-1'}>
                    <Button variant={'outline-secondary'} style={{pointerEvents: 'none'}}>
                      {item.term}
                    </Button>
                    <Button variant={'outline-secondary'} onClick={() => useArticleStore.getState().removeArticleTerm(item.term, item.type)}>
                      <X />
                    </Button>
                  </ButtonGroup>
                );
              })}
            </h5>
          </div>
        </Col>
      </Row>
      <Accordion>
        {relevantArticles.map((article, index) => (
          <Accordion.Item key={index} eventKey={index.toString()}>
            <Accordion.Header>
              <div>
                <div className={'mb-2'}>{article.title}</div>
                <div className={'d-flex align-items-center small text-secondary'}>
                  {article.authors && article.authors.join(', ')}
                  <div className={'vr mx-2'} />
                  {article.country}
                  <div className={'vr mx-2'} />
                  {article.pm_year}
                </div>
              </div>
            </Accordion.Header>
            <Accordion.Body>
              <div>{highlightWords(article.abstract, terms)}</div>
              <Row className={'mt-3'}>
                <Col className={'d-flex flex-column'}>
                  <div className={'d-flex align-items-center mb-2'}>
                    <span className={'text-secondary'}>Published in</span>
                    <div className={'vr mx-2'} />
                    {article.venue_title}
                  </div>
                  <div className={'d-flex align-items-center'}>
                    <span className={'text-secondary'}>Publication year</span>
                    <div className={'vr mx-2'} />
                    {article.venue_year}
                  </div>
                  {article.key_words.length !== 0 && (
                    <div className={'d-flex align-items-center flex-wrap mt-2'}>
                      {article.key_words.map((word, index) => (
                        <Badge key={index} bg={'secondary'} className={'me-2 my-1'}>
                          {word}
                        </Badge>
                      ))}
                    </div>
                  )}
                </Col>
                <Col className={'d-flex justify-content-end align-items-center'}>
                  <Button href={article.url} target={'_blank'} variant={'outline-primary'}>
                    <div className={'d-flex align-items-center'}>
                      <LinkSimple weight={'light'} />
                      <div className={'vr mx-2'} />
                      Read on PubMed
                    </div>
                  </Button>
                </Col>
              </Row>
            </Accordion.Body>
          </Accordion.Item>
        ))}
      </Accordion>
    </div>
  );
};

export default RelevantArticles;
